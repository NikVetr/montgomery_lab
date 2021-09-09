library(RSelenium)
library(rvest)

driver <- rsDriver(browser=c("firefox"), port = 1234L)
remote_driver <- driver[["client"]]
remote_driver$open()

url <- paste0("https://reactome.org/PathwayBrowser/")
remote_driver$navigate(url)

#wait for page to load
Sys.sleep(3)

#snag html
html <- remote_driver$getPageSource()[[1]]
html <- minimal_html(html)

webElems <- remote_driver$findElements(using = "css selector", "img[src='https://reactome.org/PathwayBrowser/Browser/clear.cache.gif']") 
cats <- as.list(html_text(html_elements(html, css = ".GOLFMBHCOSB")))
names(webElems) <- names(cats) <- cats


recursive_clicker <- function(existing_elements, element_names){
  
  for(ei in names(existing_elements)){
    
    #click the button
    webElem <- existing_elements[[ei]]
    webElem$clickElement()
    Sys.sleep(1)
    
    #get all info
    html <- remote_driver$getPageSource()[[1]]
    html <- minimal_html(html)
    new_elements <- remote_driver$findElements(using = "css selector", "img[src='https://reactome.org/PathwayBrowser/Browser/clear.cache.gif']") 
    cats <- as.list(html_text(html_elements(html, css = ".GOLFMBHCOSB")))
    names(cats) <- cats
    
    #check if we have reached the bottom level
    if(length(old_elements) == length(new_elements)){
      return(setdiff(cats, element_names))
    }
    
    names(new_elements) <- cats
    
    #subset info to just the novel elements
    new_elements <- new_elements[!sapply(names(new_elements), function(x) x %in% names(existing_elements))]
    new_element_names <- names(new_elements)
    
    print(new_element_names)
    
    #call yourself
    element_names[[ei]][[length(element_names[[ei]]) + 1]] <- recursive_clicker(existing_elements = new_elements, element_names = new_element_names)
    
  }
  
  return(element_names)
  
}


#doesn't actually work just yet lol
recursive_clicker(existing_elements = webElems, element_names = cats)

#can actually just load the pathways in directly! doh
pathways <- fread("~/data/smontgom/ReactomePathways.txt")
relation <- fread("~/data/smontgom/ReactomePathwaysRelation.txt")
str(pathways)
