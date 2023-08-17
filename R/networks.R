#' @importFrom data.table fread
#' @export retrieve_tflinks
retrieve_tflinks <- function(experimental_data=c("SS", 
                                     "LS", 
                                     "All"),
                             species=c("Homo_sapiens",
                                       "Mus_musculus", 
                                       "Rattus_norvegicus",
                                       "Drosophila_melanogaster",
                                       "Caenorhabditis_elegans", 
                                       "Saccharomyces_cerevisiae",
                                       "Danio_rerio"),
                             version="v1.0",
                             force=FALSE){
  
  print_msg("Retrieving interactions from TFLink gateway.", msg_type = "DEBUG")
  
  experimental_data <- match.arg(experimental_data)
  species <- match.arg(species)

  ext <- ifelse(experimental_data %in% c("LS", "All"), ".gz", "")
  
  tflink_basename <- paste0(species, 
                            "_interactions_", 
                            experimental_data,
                            "_simpleFormat_v1.0.tsv", 
                            ext)
  old_path <- getwd()
  dir_path <- file.path(path.expand('~'), ".scigenex", "tflink")
  dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
  setwd(dir_path)
  
  if(!file.exists(tflink_basename) | force){
    if(force)
      unlink(tflink_basename)
    download.file(paste0("https://cdn.netbiol.org/tflink/download_files/TFLink_", tflink_basename),
                  destfile=tflink_basename)

  }else{
    print_msg("Using tflink local file.", msg_type = "INFO")
  }

  interaction_full <- data.table::fread(tflink_basename, 
                                        header=TRUE,
                                        sep="\t")

  setwd(old_path) 

  return(interaction_full)
}

#' @returns a ggplot diagram.
#' @examples 
#' # load a dataset
#' load_example_dataset('7871581/files/pbmc3k_medium_clusters')
#' create_network(pbmc3k_medium_clusters)
#' @importFrom igraph graph_from_data_frame
#' @importFrom   ggraph ggraph geom_edge_link geom_node_point
#' @export create_network
create_network <- function(cluster_set=NULL,
                           experimental_data=c("SS", 
                                               "LS", 
                                               "All"),
                           species=c("Homo_sapiens",
                                     "Mus_musculus", 
                                     "Rattus_norvegicus",
                                     "Drosophila_melanogaster",
                                     "Caenorhabditis_elegans", 
                                     "Saccharomyces_cerevisiae",
                                     "Danio_rerio"),
                           version="v1.0"){
  
  check_format_cluster_set(cluster_set)
  
  cluster_set <- cluster_set[1, ]
  cluster_set_clust <- cluster_set@gene_clusters[[1]]

  print_msg("The first cluster will be analyzed by default.", msg_type = "INFO")
  
  tf_link <- retrieve_tflinks(experimental_data=experimental_data,
                              species=species,
                              version=version)
  
  test <- tf_link$Name.TF %in% cluster_set_clust & tf_link$Name.Target %in% cluster_set_clust 
  tf_link <- tf_link[test , c("Name.TF", "Name.Target", "Detection.method")]
  return(tf_link)
  graph <- igraph::graph_from_data_frame(tf_link)
  ggraph::ggraph(graph) + 
    ggraph::geom_edge_link(aes(colour = factor(year))) + 
    ggraph::geom_node_point()
}




display_network_shiny <-function(){
  
  myTag <- htmltools::tags$div(
    class = "divclass", 
    id = "first",
    htmltools::tags$h1("My first child!"),
    span(class = "child", id = "baby", "Crying")
  )
  
  myTag$children[[2]]$attribs$class <- "adult"
  
  tag_list_1 <- htmltools::tagList(
    htmltools::tags$p("Some text"),
    htmltools::tags$div("Content")
  )
  
  customTag <- htmltools::tag(
    "test", 
    list(class = "test", htmltools::tags$p("Custom Tag"))
  )
  
  
  ui <- shiny::fluidPage(
    myTag,
    HTML("<div id='chart' ></div>"),
    tag_list_1,
    htmltools::withTags(
      htmltools::tags$nav(htmltools::tags$div(), 
                     htmltools::tags$ul(htmltools::tags$li(), 
                                        htmltools::tags$li()))
    ),
    htmltools::tags$nav("This is the navigation"),
    htmltools::tags$style("p { color: red;}"),
    htmltools::tags$h1("This is a title"),
    htmltools::tags$script("
                                   var w = 960,
            h = 500,
            fill = d3.scale.category10(),
            nodes = d3.range(9).map(Object);
        
        var groups = d3.nest().key(function(d) { return d & 3; }).entries(nodes);
        
        var groupPath = function(d) {
            var fakePoints = [];
            if (d.values.length == 2)
            {
                //[dx, dy] is the direction vector of the line
                var dx = d.values[1].x - d.values[0].x;
                var dy = d.values[1].y - d.values[0].y;
        
                //scale it to something very small
                dx *= 0.00001; dy *= 0.00001;
        
                //orthogonal directions to a 2D vector [dx, dy] are [dy, -dx] and [-dy, dx]
                //take the midpoint [mx, my] of the line and translate it in both directions
                var mx = (d.values[0].x + d.values[1].x) * 0.5;
                var my = (d.values[0].y + d.values[1].y) * 0.5;
                fakePoints = [ [mx + dy, my - dx],
                              [mx - dy, my + dx]];
                //the two additional points will be sufficient for the convex hull algorithm
            }
        
            //do not forget to append the fakePoints to the input data
            return 'M' + 
                d3.geom.hull(d.values.map(function(i) { return [i.x, i.y]; })
                             .concat(fakePoints))
                .join('L') 
                + 'Z';
        }
        
        var groupFill = function(d, i) { return fill(i & 3); };
        
        var vis = d3.select('#chart').append('svg')
        .attr('width', w)
        .attr('height', h);
        
        var force = d3.layout.force()
        .nodes(nodes)
        .links([])
        .size([w, h])
        .start();
        
        var node = vis.selectAll('circle.node')
        .data(nodes)
        .enter().append('circle')
        .attr('class', 'node')
        .attr('cx', function(d) { return d.x; })
        .attr('cy', function(d) { return d.y; })
        .attr('r', 8)
        .style('fill', function(d, i) { return fill(i & 3); })
        .style('stroke', function(d, i) { return d3.rgb(fill(i & 3)).darker(2); })
        .style('stroke-width', 1.5)
        .call(force.drag);
        
        vis.style('opacity', 1e-6)
        .transition()
        .duration(1000)
        .style('opacity', 1);
        
        force.on('tick', function(e) {
          
          // Push different nodes in different directions for clustering.
          var k = 6 * e.alpha;
          nodes.forEach(function(o, i) {
            o.x += i & 2 ? k : -k;
            o.y += i & 1 ? k : -k;
          });
          
          node.attr('cx', function(d) { return d.x; })
          .attr('cy', function(d) { return d.y; });
          
          vis.selectAll('path')
          .data(groups)
          .attr('d', groupPath)
          .enter().insert('path', 'circle')
          .style('fill', groupFill)
          .style('stroke', groupFill)
          .style('stroke-width', 40)
          .style('stroke-linejoin', 'round')
          .style('opacity', .2)
          .attr('d', groupPath);
        });
        
        d3.select('body').on('click', function() {
          nodes.forEach(function(o, i) {
            o.x += (Math.random() - .5) * 40;
            o.y += (Math.random() - .5) * 40;
          });
          force.resume();
        });
        "
                
    ),
    htmltools::p(id = "hello", onclick="changeColor('green')", "Hello World")
    )
  
  server <- function(input, output, session) {}
  
  shiny::shinyApp(ui, server)
  
}
