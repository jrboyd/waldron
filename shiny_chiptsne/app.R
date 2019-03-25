#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
setwd("~/R/waldron/shiny_chiptsne")
source("app_setup.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("ChIP-tsne development"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            radioButtons("globalViewType", label = "Global View", 
                         choices = c(GLOBAL_VIEW_POINTS, 
                                     GLOBAL_VIEW_PROFILES_FAST, 
                                     GLOBAL_VIEW_PROFILES_SLOW, 
                                     GLOBAL_VIEW_DENSITY)),
            
            sliderInput("bins",
                        "Number of bins:",
                        min = 2,
                        max = 30,
                        value = 16),
            sliderInput("xrng", "x-range", -.5, .5, value = c(-.5, .5), dragRange = TRUE),
            sliderInput("yrng", "y-range", -.5, .5, value = c(-.5, .5), dragRange = TRUE),
            actionButton("doZoom", label = "Zoom"),
            selectInput("selCells", "Select Cells", choices = UI_CELLS, selected = UI_CELLS[1:4], multiple = TRUE),
            selectInput("selGenes", "Select Genes", choices = UI_GENES, selected = "RUNX1")
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            fluidRow(
                plotOutput("globalPlot", width = "400px", height = "400px") ,
                plotOutput("zoomPlot", width = "400px", height = "400px")
            ),
            fluidRow(
                plotOutput("genePlot", width = "400px", height = "400px"),
                plotOutput("profilePlot", width = "400px", height = "400px")
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$globalPlot <- renderPlot({
        typ = input$globalViewType
        if(typ == GLOBAL_VIEW_POINTS){
            p = ggplot(tsne_tp, aes(x = tx, y = ty)) + 
                geom_point() 
        }else if(typ == GLOBAL_VIEW_DENSITY){
            p = ggplot(tsne_res, aes(x = tx, y = ty)) + 
                geom_density2d()
        }else if(typ == GLOBAL_VIEW_PROFILES_FAST){
            p = ggplot(glyph_df, aes(gx, gy, group = paste(gid, mark), color = mark)) + 
                geom_path() +
                scale_color_manual(values =  c("input" = "blue",
                                               "H3K4me3" = "forestgreen",
                                               "H3K27me3" = "red")) 
        }else if(typ == GLOBAL_VIEW_PROFILES_SLOW){
            p = p_basic
        }
        
        if(any(input$xrng != c(-.5, .5)) | any(input$yrng != c(-.5, .5))){
            p = p + annotate("rect", 
                             xmin = min(input$xrng), xmax = max(input$xrng),
                             ymin = min(input$yrng), ymax = max(input$yrng),
                             fill = "#00FF0055", color = "black")           
        }
        p + 
            coord_fixed() + theme_classic() + labs(x = "", y = "")    
        
    })
    
    output$zoomPlot <- renderPlot({
        zimg_res = make_tsne_img(tsne_input$bw_dt[cell %in%  input$selCells], 
                                 tsne_res, n_points = input$bins, 
                                 xrng = zoom_xrng(), yrng = zoom_yrng())
        p = make_img_plots(zimg_res, min_size = .3, qcell = input$selCells,
                           xrng = zoom_xrng(), 
                           yrng = zoom_yrng(), 
                           as_facet = TRUE ) +
            theme_classic()
        p
    })
    
    output$imgPlot <- renderPlot({
        p_basic + 
            coord_fixed() + theme_classic() + labs(x = "", y = "")
    })
    
    output$genePlot <- renderPlot({
        plot_velocity_arrows_selected(tsne_res, 
                                      tsne_input$query_gr, 
                                      input$selCells, 
                                      tss_ids = input$selGenes) + 
            coord_fixed() + theme_classic() + labs(x = "", y = "")
    })
    
    output$profilePlot <- renderPlot({
        plot_profiles_selected(tsne_input$bw_dt, 
                               tsne_input$query_gr, 
                               input$selCells, 
                               tss_ids = input$selGenes) + 
            theme_classic() + labs(x = "", y = "")
    })
    
    zoom_xrng = reactiveVal(c(-.5, .5))
    zoom_yrng = reactiveVal(c(-.5, .5))
    
    observeEvent(input$doZoom, {
        req(input$xrng)
        req(input$yrng)
        if(any(input$xrng != zoom_xrng())){
            zoom_xrng(input$xrng)
        }
        if(any(input$yrng != zoom_yrng())){
            zoom_yrng(input$yrng)
        }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

