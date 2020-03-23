#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

library(tidyverse)
library(lubridate)
library(rstan)
options(mc.cores = parallel::detectCores())
library(bayesplot)


confirmed_sheet<-read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv") %>%
    select(-Lat, -Long)
deaths_sheet <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv")%>%
    select(-Lat, -Long)
recovered_sheet <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Recovered.csv")%>%
    select(-Lat, -Long) 



confirmed_long <- gather(confirmed_sheet, -Province.State , -Country.Region , key="date", value = "cases")%>%
    separate(date,c("x","date.longer"),1,remove=TRUE) %>% 
    separate(date.longer,c("month","day","year"),remove=TRUE) %>%
    separate(Province.State,c("location","State"),sep=", ",remove=FALSE) %>% #for US data before 3/10/2020 Province.State value includes sub-state location data. Split this out so that we can recover state level data later.
    mutate(location = as.character(location)) %>%
    mutate(State = as.character(State)) %>%
    mutate(year=as.character(as.numeric(year)+2000)) %>% #data was in a format with just the last two digits for year
    unite(comb.date, c(month,day,year) , sep=".")%>%  
    mutate(date = parse_date(comb.date , "%m%.%d%.%Y"))%>%
    select(-comb.date , -x) %>%
    mutate(delta.days=period_to_seconds(days(ymd(date) - ymd(20200122)))/(60*60*24)) %>% #calculate days since data began being collected
    as_tibble()  


maxdays=max(confirmed_long$delta.days)




USTotals <- filter(confirmed_long,Country.Region=="US" ) %>%  
    group_by(date, delta.days, Country.Region) %>% summarise(mean=mean(cases), count=n() ) %>%
    mutate(total.cases = mean*count) %>% 
    select(-mean, -count)  %>%
    #  mutate(Country.Region = "USA") %>%
    rename(cases=total.cases) 

ChinaTotals <- filter(confirmed_long,Country.Region=="China" ) %>%  
    group_by(date, delta.days, Country.Region) %>% summarise(mean=mean(cases), count=n() ) %>%
    mutate(total.cases = mean*count) %>% 
    select(-mean, -count)  %>%
    mutate(Country.Region = "China") %>%
    rename(cases=total.cases) 


CanadianTotals <- filter(confirmed_long,Country.Region=="Canada" ) %>%  
    group_by(date, delta.days, Country.Region) %>% summarise(mean=mean(cases), count=n() ) %>%
    mutate(total.cases = mean*count) %>% 
    select(-mean, -count)  %>%
    rename(cases=total.cases) 

confirmed_long2 <- bind_rows(filter(confirmed_long,Country.Region!="Canada",Country.Region!="China",Country.Region!="US"),USTotals,ChinaTotals,CanadianTotals) 



#make list of countries by most cases
country_list <- confirmed_long2 %>% filter(delta.days==maxdays) %>% arrange(desc(cases)) %>% select(Country.Region)
#countries in alphabetical order but only with more than 100 cases
country_list_alpha <- confirmed_long2 %>% filter(delta.days==maxdays,cases>100) %>% arrange(Country.Region) %>% select(Country.Region)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Compare two regions"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("country1", "Country 1:",
                        country_list),
            selectInput("country2", "Country 2:",
                        country_list_alpha),
            
            sliderInput("start1",
                        "Start day 1:",
                        min = 1,
                        max = 60,
                        value = 25),
            sliderInput("end1",
                        "End day 1:",
                        min = 1,
                        max = 60,
                        value = 60),
            sliderInput("start2",
                        "Start day 2:",
                        min = 1,
                        max = 60,
                        value = 25),
            sliderInput("end2",
                        "End day 2:",
                        min = 1,
                        max = 60,
                        value = 60),
            actionButton("do", "Run Models")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("country1"),plotOutput("country2"),plotOutput("compare")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    observeEvent(input$do, {
        tocomp = tibble(country = c(input$country1,input$country2) , startday = c(input$start1,input$start2) , 
                        endday = c(input$end1,input$end2)) %>%
        mutate(n=endday-startday)
    
    
    mydata.sub1 <- filter(confirmed_long2, delta.days>=tocomp$startday[1],delta.days<=tocomp$endday[1],Country.Region == tocomp$country[1]) %>% mutate(dataset=1)   
    mydata.sub2 <- filter(confirmed_long2, delta.days>=tocomp$startday[2],delta.days<=tocomp$endday[2],Country.Region == tocomp$country[2])   %>% mutate(dataset=2)
    
    mindays=c(min(mydata.sub1$delta.days),min(mydata.sub2$delta.days))
    
    mydata.sub1<-mydata.sub1 %>% mutate(days=delta.days-mindays[1]+1)
    mydata.sub2<-mydata.sub2 %>% mutate(days=delta.days-mindays[2]+1)
    
    maxdays=c(max(mydata.sub1$days),max(mydata.sub2$days))
    
    mydata.sub.total <-bind_rows(mydata.sub1,mydata.sub2) %>% group_by(dataset) %>% arrange(days, .by_group = TRUE) %>%
        rename(total_cases=cases) 
    
    total_datapoints=nrow(mydata.sub.total)
    
    stan_data <- c(mydata.sub.total[c("days","total_cases","dataset")],tocomp[c("country")],list(maxdays=maxdays,total_datapoints=total_datapoints , pstrength =0.5)) 
    
    
    fit_poissCompare <- stan(file = 'model_compare2.stan', 
                             data =stan_data, chains = 8,iter = 5000, seed = 2131231 )
    
    
    posterior <- as.matrix(fit_poissCompare)
    
    output$compare <-  renderPlot({mcmc_areas(posterior,
               pars = c("lambda_diff"),
               prob = 0.99) })
    })
    

    output$country1 <- renderPlot({
        ggplot(data=filter(confirmed_long2,  Country.Region==input$country1  ,delta.days>input$start1  ,delta.days<input$end1 ,cases>20) , aes(x=delta.days,y=log(cases,base=10))) +
            geom_point()+
            geom_smooth(method = "lm", formula = y ~ x) + 
            scale_y_continuous( limits=c(1,6),breaks=c(1,2,3,4,5), labels=c(10,100,1000,10000,100000))+
            scale_x_continuous( limits=c(20,56))+
            labs( x="days since Jan 22" , y="cases", title=input$country1)
    })
    
    output$country2 <- renderPlot({
        ggplot(data=filter(confirmed_long2,  Country.Region==input$country2 ,delta.days>input$start2  ,delta.days<input$end2  ,cases>20) , aes(x=delta.days,y=log(cases,base=10) )) +
            geom_point()+
            geom_smooth(method = "lm", formula = y ~ x) + 
            scale_y_continuous( limits=c(1,6),breaks=c(1,2,3,4,5), labels=c(10,100,1000,10000,100000))+
            scale_x_continuous( limits=c(20,56))+
            labs( x="days since Jan 22" , y="cases", title=input$country2)
    })
    
    
    

        
}

# Run the application 
shinyApp(ui = ui, server = server)
