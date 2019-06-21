library(shiny)

Sweave('App_mapGenesOnSingleCellClusters.Rnw')

## UI
ui <- fluidPage(
    titlePanel(h3('Single Cell Analysis (In vitro differentiation of hESCs into pancreatic endocrine cells, based on: Sharon et al., 2019, Cell Reports 27, 2281–2291 May 21, 2019 https://doi.org/10.1016/j.celrep.2019.04.083'), windowTitle='SingleCell_MEP'),
    sidebarLayout(
        sidebarPanel(
            conditionalPanel(
                condition = 'input.thetabs == 1',
                textInput(inputId='symbol', 'Symbol of Gene to map on clusters', value='VIM'),
                checkboxInput('log', 'Log transform?', value=TRUE, width=NULL),
                verbatimTextOutput('info')),
            conditionalPanel(
                condition = 'input.thetabs == 2',
                selectInput('eprog', 'Choose an EP' , c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',
                                                        '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25',
                                                        '26','27', '28', '29', '30', '31', '32', '33', '34', '35', '36')),
                sliderInput('psize', 'Point Size', min=0.5, max=4, value=3)),
            conditionalPanel(
                condition = 'input.thetabs == 3',
                checkboxGroupInput('chckeps', 'Choose one or more EPs', choices= c(
                                                                            EP1 = 'EP_1', EP2 = 'EP_2', EP3 = 'EP_3', EP4 = 'EP_4', EP5 = 'EP_5',
                                                                            EP6 = 'EP_6', EP7 = 'EP_7', EP8 = 'EP_8', EP9 = 'EP_9', EP10 = 'EP_10',
                                                                            EP11 = 'EP_11', EP12 = 'EP_12', EP13 = 'EP_13', EP14 = 'EP_14',
                                                                            EP15 = 'EP_15', EP16 = 'EP_16', EP17 = 'EP_17', EP18 = 'EP_18',
                                                                            EP19 = 'EP_19', EP20 = 'EP_20', EP21 = 'EP_21', EP22 = 'EP_22',
                                                                            EP23 = 'EP_23', EP24 = 'EP_24', EP25 = 'EP_25', EP26 = 'EP_26',
                                                                            EP27 = 'EP_27', EP28 = 'EP_28', EP29 = 'EP_29', EP30 = 'EP_30',
                                                                            EP31 = 'EP_31', EP32 = 'EP_32', EP33 = 'EP_33', EP34 = 'EP_34',
                                                                            EP35 = 'EP_35', EP36 = 'EP_36'), selected = c('EP_6'), inline=TRUE),
                sliderInput('mpsize', 'Point Size', min=0.5, max=4, value=1.5)),
            conditionalPanel(
                condition = 'input.thetabs == 4',
                checkboxGroupInput('chck4', 'Choose one or more EPs', choices= c(
                                                                            EP1 = 'EP_1', EP2 = 'EP_2', EP3 = 'EP_3', EP4 = 'EP_4', EP5 = 'EP_5',
                                                                            EP6 = 'EP_6', EP7 = 'EP_7', EP8 = 'EP_8', EP9 = 'EP_9', EP10 = 'EP_10',
                                                                            EP11 = 'EP_11', EP12 = 'EP_12', EP13 = 'EP_13', EP14 = 'EP_14',
                                                                            EP15 = 'EP_15', EP16 = 'EP_16', EP17 = 'EP_17', EP18 = 'EP_18',
                                                                            EP19 = 'EP_19', EP20 = 'EP_20', EP21 = 'EP_21', EP22 = 'EP_22',
                                                                            EP23 = 'EP_23', EP24 = 'EP_24', EP25 = 'EP_25', EP26 = 'EP_26',
                                                                            EP27 = 'EP_27', EP28 = 'EP_28', EP29 = 'EP_29', EP30 = 'EP_30',
                                                                            EP31 = 'EP_31', EP32 = 'EP_32', EP33 = 'EP_33', EP34 = 'EP_34',
                                                                            EP35 = 'EP_35', EP36 = 'EP_36'), selected = c('EP_6'), inline=TRUE),
                #textInput(inputId='epnum', 'Integers separated by commas', 6),
                textInput(inputId='sym4', 'Symbole of Gene to Boxplot', 'VIM'),
                actionButton('abutton', 'GO'))
        ),
        mainPanel(
            tabsetPanel(id='thetabs', selected=1, type='tabs',
                        tabPanel('Simlr', value=1, fluidRow(column(12, offset=1, plotOutput('mapping', width='700px', height='600px')),
                                                            column(12, offset=1, plotOutput('stages', width='700px', height='600px')))),
                        tabPanel('EP', value=2, fluidRow(column(12, offset=1, plotOutput('ep', width='700px', height='600px')))),
                        tabPanel('MultEP', value=3, fluidRow(column(12, offset=1, plotOutput('mep', width='700px', height='600px')))),
                        tabPanel('Boxplot', value=4, fluidRow(column(12, offset=1, plotOutput('bxp', width='700px', height='600px'))))
                        )
        )
    )
)

## Server
server <- function(input, output) {
    output$mapping <- renderPlot({
        output$info <- renderText({'Only One Symbol'})

        symreact <- reactive({
            validate(
                need(any(grepl(toupper(sub('\\.', '\\\\.', input$symbol)), ls(sym2name))) & input$symbol != '',
                     "This symbol is not recognized or the gene is not expressed in any of the cells! Please, choose another gene.")
            )
            toupper(input$symbol)
        })
        
        #sym <- toupper(symreact())
        if (input$log) {
            y.df <- mergeScoresAndGeneExpr(10, lset, symreact())
            thistitle <- 'FPKM\n(log2)'
        } else {
            #y.df <- mergeScoresAndGeneExpr(10, nset, symreact())
            y.df <- mergeScoresAndGeneExpr(10, nlset, symreact())
            thistitle <- 'FPKM'
        }

        name <- mget(symreact(), sym2name, ifnotfound=NA)
        name <- as.character(unlist(name))
        output$info <- renderText({name})
        
        p <- ggplotGenesOnCells(y.df, symreact(), 10, thistitle)
        plot(p)
    })
    
    output$stages <- renderPlot({
        pstage <- ggplotStagesOnSIMLRclusters(nsetH10c.simlr, nset)
        plot(pstage)
    })

    epreact <- reactive({
        paste('EP', input$eprog, sep='_')
    })

    output$ep <- renderPlot({
        y.df <- mergeScoresAndGeneExpr_vEP(nset, nsetH10c.simlr, acep.df, input$eprog)
        p <- ggplotGenesOnCells_vEP(y.df, epreact(), 10, input$psize)
        print(p)
    })

### multEP tab ###
    output$mep <- renderPlot({
        pmap <- ggplotMultEPonSIMLRclusters(nset, nsetH10c.simlr, acep.df, input$mpsize, input$chckeps)
        print(pmap)
    })

    ensid <- eventReactive(input$abutton, {
        sl <- input$sym4
        csl <- toupper(sub('\\.', '\\\\.', input$sym4))
        as.character(unlist(mget(csl, sym2ens, ifnotfound=NA)))
    })
        
    output$bxp <- renderPlot({
        sel <- input$chck4
        df <- selMaxEP_vExt(acep.df, sel)
        
        bx <- ggboxplotOverCells(df, lset, ensid())
        print(bx)
    })
}

shinyApp(ui=ui, server=server)
    
