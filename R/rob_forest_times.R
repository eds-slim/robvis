#' @title Append a risk-of-bias traffic-light plot to a forest plot
#'
#' @description A wrapper for metafor::forest function, which adds a risk of
#'   bias traffic-light plot to the right-hand side of the forest plot. The
#'   heavy lifting for this function is done by metafor. Note that if not
#'   specified as additional arguments, this functions sets the header argument
#'   of metafor::forest() to TRUE.
#'
#' @param res Output from metafor meta-analysis function
#' @param rob_tool The risk-of-bias assessment tool used to perform the
#'   assessments
#' @param rob_me Optional value defining the result of the Risk-Of-Bias due to
#'   Missing Evidence (ROB-ME) assessment for this synthesis. By default (rob_me
#'   = NULL), this is omitted from the plot.
#' @param rob_levels Vector otoof judgments [e.g. c("Low","Some
#'   concerns","High","Critical")] that controls the ordering of subgroups
#'   within the plot
#' @param title Text to use for plot title
#' @param rob_legend Logical specifying whether a legend for the risk-of-bias
#'   plot should be shown. Default is TRUE.
#' @param rob_legend_cex Expansion factor for the risk-of-bias legend
#' @param ... Additional arguments to be passed to the metafor::forest()
#'   function
#'
#' @family main
#'
#' @export

rob_forest_times <-
  function(res
           , rob_me = NULL
           , rob_levels = NULL
           , title = NULL
           , rob_legend = TRUE
           , rob_legend_cex = 0.9
           , group.var = 'overall'
           , layout = c(1,1,1)
           , es = 'OTG/ min'
           , ...
           ) {

    # mod by ES
    
    # Check that res is of class RMA
    if (!("rma" %in% class(res))) {
      stop("Result objects need to be of class \"meta\" - output from metafor package functions")
    }


    colnames(res$data$df) <- stringr::str_to_lower(colnames(res$data$df))


    dat <- res$data$df |> 
      dplyr::mutate(overall = factor(overall, levels = rob_levels)) |> 
      dplyr::arrange(!!sym(group.var), desc(year))
       


    # Get maximum domain

    max_domain_column <- dat |> 
      dplyr::select(dplyr::matches("^d.$")) |> 
      colnames() |> 
      gsub("d", "", x = _) |> 
      as.numeric() |> 
      max() + 2


    # Use this to define the gaps between different groups
    # Will be important when adding argument to prevent subgroup analyses
    offset_n <- 3

    dat_rob_vec <- dat |> 
      dplyr::mutate(row_n = 1:dplyr::n()) |> 
      dplyr::group_by(!!sym(group.var)) |> 
      dplyr::summarise(n = dplyr::n(), max = max(row_n), min = min(row_n)) %>% 
      dplyr::mutate(offset = seq(1, length(unique(.[[group.var]]))*offset_n, by = offset_n)) |> 
      dplyr::mutate(min = min+offset, max = max+offset, heading = max+1, stats = min-1.25) |> 
      dplyr::mutate(min = ifelse(n==1, min-1, min)
                    , max = ifelse(n==1, max-1, max)
                    , heading = ifelse(n==1, heading-1, heading))

    if (length(unique(dat[[group.var]]))==1) {
      dat_rob_vec <- dat_rob_vec |> 
        dplyr::mutate(dplyr::across(c(min, max, heading), ~. -1))
    }
    
    #res <- stats::update(res, data = dat) ## no update fpr metamedian analysis

    
    rows <- c()

    for (i in 1:nrow(dat_rob_vec)) {

      rows <-c(rows, dat_rob_vec$min[i]:dat_rob_vec$max[i])

    }
    rows <- rows + 1
    
    print(rows)

    arg <- list(...)

    dd <- (arg$at[length(arg$at)] - arg$at[1])/layout[2]
    x_min <- arg$at[1] - dd*layout[1]
    arg$ilab.xpos <- arg$at[1] - c(1/3+1/5, 1/3, 1/5, 0)*dd*layout[1]
   
  

    x_max = arg$at[length(arg$at)] + dd*layout[3]/3
    textpos <- c(x_min, x_max)
    y_max <- max(rows)+4

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    # Deal with adding rob data

    dat <- dat |> 
      dplyr::mutate(dplyr::across(dplyr::matches("^d.$|overall"), clean_data))

    ddd <- 2/3*dd*layout[3]/(max_domain_column)
    x_pos <- seq(x_max + ddd, by = ddd, length.out = max_domain_column - 2)

    x_overall_pos <- max(x_pos) + ddd

    # Convenience vector, specifying x-axis positions for all risk of bias columns
    header_row <- c(x_pos, x_overall_pos)

    legend_pos <- x_max+(max(header_row)-min(header_row))/2

    # New right-hand x-axis limit
    new_x_lim <- x_overall_pos

    rob_colours <- get_colour('ROBINS-I', "colourblind")


      judgements<-   c("Critical risk of bias"
                       , "Serious risk of bias"
                       , "Moderate risk of bias"
                       , "Low risk of bias"
                       , "No information"
                       )
      cols <- c(
        c = rob_colours$critical_colour
        , s = rob_colours$high_colour
        , m = rob_colours$concerns_colour
        , l = rob_colours$low_colour
        , n = rob_colours$ni_colour
        , x = rob_colours$na_colour
      )

      syms <- c(c = "!"
                , s = "×"
                , m = "?"
                , l = "+"
                , n = ""
                , x = ""
                )


      shapes <- c(c = 19
                  , s = 19
                  , m = 19
                  , l = 19
                  , n = 19
                  , x = 19
                  )

    rob_psize = 2.5
    tsize <- rob_psize * 0.3

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    # Make forest plot

    if (is.null(arg$header)) {
      arg$header = "Author(s) and Year"
    }

    if (is.null(arg$addpred)) {
      arg$addpred = FALSE
    }

    if (is.null(arg$mlab)) {
      arg$mlab = mlabfun("RE Model for all studies", res)
    }
    

    arg$x <- res
    arg$xlim <- c(x_min, new_x_lim)
    arg$ylim=c(-2, y_max)
    arg$rows <- rows
    arg$textpos <- textpos
    arg$col <- 'darkblue'
    arg$ilab <- cbind(dat$n.g1, dat$es.g1, dat$n.g2, dat$es.g2)
    arg$ilab.lab <- c('n', es, 'n', es)
    arg$ilab.pos <- 2

    ### set up forest plot (with 2x2 table counts added; the 'rows' argument is
    ### used to specify in which rows the outcomes will be plotted)
    # metafor::forest(res, xlim=c(x_min, new_x_lim), atransf=exp,
    #        cex=1.2, ylim=c(-1.5, y_max), rows=rows, textpos=textpos,
    #        mlab=mlab, addpred = addpred)
    
    par(mar = c(5, 4, 4, 2) + 0.1, xpd = NA)
    f <- do.call(metafor::forest, arg)

    graphics::text(mean(f$ilab.xpos[2]), y_max, 'Direct transport', font=2, adj = 1)
    graphics::text(mean(f$ilab.xpos[4]), y_max, 'Control', font=2, adj = 1)
    
    
    segments(x0 = f$ilab.xpos[1], y0 = y_max - .5, x1 = f$ilab.xpos[2], y1 = y_max - .5 )
    segments(f$ilab.xpos[3], y_max - .5, f$ilab.xpos[4], y_max - .5)
    
    
    ### set font expansion factor (as in forest() above) and use a bold font

    op <- graphics::par(font=2)

    ### switch to italic font
    graphics::par(font = 3)

    ### add text for the subgroups
    for (i in 1:nrow(dat_rob_vec)) {

      graphics::text(x_min, dat_rob_vec$heading[i] + 1, pos = 4, dat_rob_vec[[group.var]][i], cex = 1.2)
    }

    ### set par back to the original settings
    graphics::par(op)

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    # Add risk of bias data

    headers <- c(paste0("D",seq_len(max_domain_column-2)),"O")

    graphics::par(font = 2)
    # Need to add handling of top here
    graphics::text(mean(header_row), y_max, labels = "Risk of Bias")
    graphics::text(header_row, y_max-2 + 1, labels = headers)
    segments(x0 = header_row[1], y0 = y_max - .5, x1 = header_row[length(header_row)], y1 = y_max - .5 )
    
    graphics::par(op)

    # Plot domain points
    for (j in 1:length(x_pos)) {
      graphics::points(
        x = rep(x_pos[j], length(rows))
        , y = rows
        , pch = shapes[dat[[paste0("d", j)]]]
        , col = scales::alpha(cols[dat[[paste0("d", j)]]], 0.6)
        , cex = rob_psize
      )
      graphics::text(x_pos[j], rows, syms[dat[[paste0("d", j)]]], cex = tsize)
    }


    graphics::points(
      rep(x_overall_pos, length(rows))
      , rows
      , pch = 15
      , col = scales::alpha(cols[dat[["overall"]]], 0.6)
      , cex = rob_psize
    )
    graphics::text(x_overall_pos, rows, syms[dat[["overall"]]], cex = tsize)
    graphics::par(op)
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    # Add sub-group, summary polygons & text

    rma_flexi <- function(x) {
        rma(yi = res$yi, vi = res$vi, subset = (res$data$df[[group.var]] == x)
        )
      }

    subgroup_res <- purrr::map(unique(dat[[group.var]]), rma_flexi)

    
    if (length(unique(dat[[group.var]])) > 1) {

      ### add summary polygons for the subgroups
      for (i in 1:nrow(dat_rob_vec)) {

        if (length(subgroup_res[[i]]$slab) == 1) {
          next
        }
        
        print(subgroup_res[[i]])
        

        metafor::addpoly(
          subgroup_res[[i]]
          #, fonts = c('serif'=3, 'mono'=3)
          , row = dat_rob_vec$stats[i] + 1
          , textpos=textpos
          , col = 'lightblue'
          ,  annotate = F
          ,  mlab = mlabfun("\tRE Model for Subgroup", subgroup_res[[i]])
        )
        
        annotate_poly(subgroup_res[[i]]$b
                      , subgroup_res[[i]]$ci.lb
                      , subgroup_res[[i]]$ci.ub
                      , textpos = textpos[2]
                      , atransf = arg$transf
                      , rows = dat_rob_vec$stats[[i]] + 1
                      , font = 3
                      )

      }
    }
  
  rect(f$textpos[2], -1.5, arg$at[length(arg$at)], -0.5, col = "white", border = NA)
  annotate_poly(res$b, res$ci.lb, res$ci.ub
                , textpos = textpos[2]
                , atransf = arg$transf
                , rows = -1
                , font = 2
                )
    
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

    if (length(unique(dat[[group.var]])) > 1 && nrow(dat) > 1) {

      # Fit meta-regression model to test for subgroup differences
      subgroup_res <- rma(yi = res$yi, vi = res$vi, mods = ~ dat[[group.var]], method = "DL")

      ### add text for the test of subgroup differences
      graphics::text(x_min, -1.8, pos = 4,  bquote(
        paste(
          "Test for Subgroup Differences: ",
          Q[M],
          " = ",
          .(formatC(
            subgroup_res$QM, digits = 2, format = "f"
          )),
          ", df = ",
          .(subgroup_res$p - 1),
          ", p = ",
          .(formatC(
            subgroup_res$QMp, digits = 2, format = "f"
          ))
        )
      ))
    }
    
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    
    if(!is.null(title)){
      graphics::par(font = 2)
      graphics::text(x_min, y_max, pos=4, bquote(bold(underline(.(title)))), cex = 1.2)
      graphics::par(op)
    }


    if (rob_legend == TRUE) {

      graphics::legend(
        legend_pos
        , -1.8
        , judgements
        , pch = 15
        , xjust = 0.5
        , col = utils::head(cols, -1)
        , xpd = TRUE
        , title = parse(text = "bold(\"Judgement\")")
        , title.adj = 0.1
        , cex = rob_legend_cex
        , pt.cex = rob_legend_cex
        , y.intersp = 0.7
      )
    }


  }
environment(rob_forest_times) <- environment(rob_forest)
