library(shiny)
ui <- fluidPage(
  titlePanel("Comparing the WL linear and the QWL model"),
  sidebarLayout(
    sidebarPanel( 
      sliderInput("a1",
                  "a1:",
                  min = -5,
                  max = 5,
                  value = 0),
      sliderInput("a2",
                  "a2:",
                  min = -5,
                  max = 5,
                  value = 0),
      sliderInput("b1",
                  "b1:",
                  min = -5,
                  max = 5,
                  value = 1),
      sliderInput("b2",
                  "b2:",
                  min = -5,
                  max = 5,
                  value = 1),
      sliderInput("mu1",
                  "mu1:",
                  min = -5,
                  max = 5,
                  value = 0),
      sliderInput("mu2",
                  "mu2:",
                  min = -5,
                  max = 5,
                  value = 0),
      sliderInput("sd1",
                  "sd1:",
                  min = 0,
                  max = 5,
                  value = 1),
      sliderInput("sd2",
                  "sd2:",
                  min = 0,
                  max = 5,
                  value = 1)
    ),
    mainPanel(
      plotOutput("distPlot", width = '600px', height = '600px'),
      pre('
For linear phi, at the COU phi1 = phi2, 
so  a1+b1*x1 = a2+b2*x2.  
If x1 and x2 are jointly normal, the corresponding CDFs are: p1 = pnorm(x1, mu1, sd1) and p2 = pnorm(x2, mu2, sd2) where
x2 =  (a1+b1*x1 – a2)/b2. Thus
           
           p2 = pnorm(x2, mu2, sd2) =  pnorm( (a1+b1*x1 – a2)/b2,     mu2, sd2)
           =  pnorm( ((a1+b1*(qnorm(p1, mu1, sd1)  – a2)/b2), mu2, sd2)

This mapping from p1 to p2 gives us a comparison to the delta model.
           
Note that cor(x1,x2) has no effect on the COU.
           
From symmetry, 
changing  a1 or mu1 makes a graph similar to the delta graphs,
while changing  b1 or sd1 makes a graph sinusoidal.
           
Note that qnorm(p1, mu1, sd1)   = qnorm(p1) sd1 + mu1
Setting b1=b2 and sd1 = sd2,
p2 = pnorm( ((a1+b1*( qnorm(p1) sd1 + mu1  – a2)/b2), mu2, sd2)
	= pnorm( (a1-a2)/b2 + ( qnorm(p1) sd1 + mu1), mu2, sd2)
	= pnorm(  ( (a1-a2)/b2 + qnorm(p1) sd1 + mu1 – mu2 )/sd2 )
	= pnorm ( ( qnorm(p1) + (mu1-mu2) +(a1-a2)/b2)/sd2 )

Defining H = qnorm and Hinv = pnorm, and 
p2 = Hinv(  H(p1) + delta )
where delta = (mu1-mu2+(a1-a2)/b2)/sd2
')
    )
  )
)
# Define server logic required to draw a histogram
server <- function(input, output) {
  output$distPlot <- renderPlot({
    a1=input$a1
    a2=input$a2
    b1=input$b1
    b2=input$b2
    mu1=input$mu1
    mu2=input$mu2
    sd1=input$sd1
    sd2=input$sd2
    print(mu1)
    print(sd1)
    p1 = seq(0, 1, length=100)
    x1 = qnorm(p1, mu1, sd1)
    x2 = ( a1+b1*x1 - a2)/b2
    p2 = pnorm( x2, mu2, sd2)
    plot(p1, p2, col='green', cex=2, pch=16)
    for(delta in seq(-4,4,length=12))
      lines(p1vec<-seq(0,1,length=100),
            deltaMap(delta = delta, p1vec))
    deltaEquiv = (a1-a2)/b2 + (mu1-mu2)/sd2
    lines(p1vec<-seq(0,1,length=100),
          deltaMap(delta = deltaEquiv, p1vec),
          col='blue', lwd=3)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

