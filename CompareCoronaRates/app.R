# A little shiny app to compare the log growth estimates of two segments of coronavirus case number data 
# Written by Stephen Proulx sproul@ucsb.edu
# library(rsconnect)
# rsconnect::deployApp('/Users/proulx/coronavirusEstimator/CoronavirusGraphs_Fits/CompareCoronaRates')

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


maxdays=max(confirmed_long$delta.days)-1




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

mod <- readRDS("model.rds") 



# Define UI
ui <- fluidPage(

    # Application title
    titlePanel("Compare two regions"),

    
    sidebarLayout(
        sidebarPanel(
            selectInput("country1", "Country 1:",
                        country_list,
                        selected = "US"),
            selectInput("country2", "Country 2:",
                        country_list_alpha,
                        selected = "Italy"),
            
            sliderInput("start1",
                        "Start day 1:",
                        min = 1,
                        max = maxdays,
                        value = 38),
            sliderInput("end1",
                        "End day 1:",
                        min = 1,
                        max = maxdays,
                        value = maxdays),
            sliderInput("start2",
                        "Start day 2:",
                        min = 1,
                        max = maxdays,
                        value = 45),
            sliderInput("end2",
                        "End day 2:",
                        min = 1,
                        max = maxdays,
                        value = maxdays),
            sliderInput("pstrength",
                        "Advanced: Strength of the prior that the two periods are the same",
                        round = -2, step = 0.01,
                        min = 0,
                        max = 2,
                        value = 0.5),
            actionButton("do", "Run Models")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("country1"),plotOutput("compare")
        )
    )
)

# Define server logic 
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
    
    stan_data <- c(mydata.sub.total[c("days","total_cases","dataset")],tocomp[c("country")],list(maxdays=maxdays,total_datapoints=total_datapoints , pstrength = input$pstrength)) 
    
 
#    fit_poissCompare <- sampling(mod, 
#                    data =stan_data, chains = 8,iter = 5000, seed = 2131231 )
                             
   fit_poissCompare <- stan(file = 'model_compare2.stan', 
                             data =stan_data, chains = 8,iter = 5000, seed = 2131231 )
    
    
    tmp<-as.data.frame(rstan::extract(fit_poissCompare,pars=c("mean_lambda","lambda_diff","lambda[1]","lambda[2]"))) %>% 
        rename("lambda1"="lambda.1.","lambda2"="lambda.2.")
    
    
    tmp1<-tmp %>%select(lambda1) %>%rename(mode=lambda1) %>%mutate(parameter="lambda1") %>% mutate(dtime=log(2)/(log(1+mode)))
    tmp2<-tmp %>%select(lambda2) %>%rename(mode=lambda2) %>%mutate(parameter="lambda2") %>% mutate(dtime=log(2)/(log(1+mode)))
    tmp<-rbind(tmp1,tmp2)
    
    
    lq=0.025
    uq=0.975
    d1<-density(tmp1$dtime)
    dd1 <- with(d1,data.frame(x,y)) %>% filter(y>0.01)
    qs1=quantile(ecdf(tmp1$dtime),prob=c(lq,0.5,uq))
    d2<-density(tmp2$dtime)
    dd2 <- with(d2,data.frame(x,y)) %>% filter(y>0.01)
    qs2=quantile(ecdf(tmp2$dtime),prob=c(lq,0.5,uq))
    
    maxdd1=max(dd1$y)
    dd1<-mutate(dd1,rely=y/maxdd1)
    maxdd2=max(dd2$y)
    dd2<-mutate(dd2,rely=y/maxdd2)
    
    
    output$compare <-  renderPlot({   
        ggplot(data=dd1 , aes(x,rely))+
            geom_line(data=dd1)+
            geom_ribbon(data=filter(dd1,x>qs1[[1]] & x<qs1[[3]]),aes(ymax=rely),ymin=0,
                        fill="#1B9E77",colour=NA,alpha=0.5)+
            geom_line(data=dd2)+
            geom_ribbon(data=filter(dd2,x>qs2[[1]] & x<qs2[[3]]),aes(ymax=rely),ymin=0,
                        fill="#D95F02",colour=NA,alpha=0.5)+
            scale_y_continuous(limits=c(0,1.5),name="Posterior Density")+
            scale_x_continuous(name="Doubling Time", limits=c(2, 8)) +
            annotate("text", x=c(qs1[[2]],qs2[[2]]),
                     y=1.2, label= c(str_c(input$country1, " 1"),str_c(input$country2, " 2")),
                     size=3)+
            theme_bw()+
            theme( axis.text=element_text(family="Helvetica", size=8),text=element_text(family="Helvetica", size=12))
    })
        
    
#    posterior <- as.matrix(fit_poissCompare)
    
#    output$compare <-  renderPlot({mcmc_areas(posterior, pars = c("lambda_diff"),prob = 0.99) })

})
    
    
    
    
    output$country1 <- renderPlot({
        data1<-filter(confirmed_long2,  Country.Region==input$country1,delta.days>(input$start1-1)  ,delta.days<(input$end1+1) ,cases>20) %>% 
            mutate(country=input$country1 ,dataset=1) %>% unite("country_set",country:dataset)
        data2<-filter(confirmed_long2,  Country.Region==input$country2  ,delta.days>(input$start2-1)  ,delta.days<(input$end2+1) ,cases>20) %>% 
            mutate(country=input$country2 ,dataset=2) %>% unite("country_set",country:dataset)
        datatot<-bind_rows(data1,data2)
        
            ggplot(data=datatot , aes(x=delta.days,y=log(cases,base=10) , color=country_set,group=country_set)) +
            geom_point()+
            geom_smooth(method = "lm", formula = y ~ x) + 
            scale_y_continuous( limits=c(1,6),breaks=c(1,2,3,4,5), labels=c(10,100,1000,10000,100000))+
            scale_x_continuous( limits=c(20,maxdays))+
            labs( x="days since Jan 22" , y="cases", title="Log case numbers")
    })

    # output$country1 <- renderPlot({
    #     ggplot(data=filter(confirmed_long2,  Country.Region==input$country1  ,delta.days>input$start1  ,delta.days<input$end1 ,cases>20) , aes(x=delta.days,y=log(cases,base=10))) +
    #         geom_point()+
    #         geom_smooth(method = "lm", formula = y ~ x) + 
    #         scale_y_continuous( limits=c(1,6),breaks=c(1,2,3,4,5), labels=c(10,100,1000,10000,100000))+
    #         scale_x_continuous( limits=c(20,56))+
    #         labs( x="days since Jan 22" , y="cases", title=input$country1)
    # })
    # 
    # output$country2 <- renderPlot({
    #     ggplot(data=filter(confirmed_long2,  Country.Region==input$country2 ,delta.days>input$start2  ,delta.days<input$end2  ,cases>20) , aes(x=delta.days,y=log(cases,base=10) )) +
    #         geom_point()+
    #         geom_smooth(method = "lm", formula = y ~ x) + 
    #         scale_y_continuous( limits=c(1,6),breaks=c(1,2,3,4,5), labels=c(10,100,1000,10000,100000))+
    #         scale_x_continuous( limits=c(20,56))+
    #         labs( x="days since Jan 22" , y="cases", title=input$country2)
    # })
    # 
    
    

        
}

# Run the application 
shinyApp(ui = ui, server = server)
