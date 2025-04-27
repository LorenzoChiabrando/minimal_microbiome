generate_pnpro <- function(arc_df, pnpro_out) {
  # Step 1: Extract unique places and transitions
  places <- unique(arc_df$place)
  transitions <- unique(arc_df$transition)
  
  # Create a lookup table for commands
  command_df <- arc_df[!is.na(arc_df$command) & arc_df$command != "", c("transition", "command")]
  # Remove duplicates (if any)
  command_df <- command_df[!duplicated(command_df$transition), ]
  # Create commands lookup
  commands <- setNames(command_df$command, command_df$transition)
  
  # Step 2: Create the XML content
  xml_header <- '<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- This project file has been saved by the New GreatSPN Editor, v.100 --><project name="Generated Project" version="121">
  <gspn name="PetriNet" zoom="50">
    <nodes>'

  # Places with automatic layout
  places_xml <- ""
  y <- 3.0
  for (place in places) {
    places_xml <- paste0(places_xml, 
                         sprintf('\n      <place label-x="0.0" label-y="0.0" name="%s" x="6.0" y="%s"/>', 
                                 place, y))
    y <- y + 2.0  # Increment y position for layout
  }
  
  # Transitions with automatic layout
  transitions_xml <- ""
  y <- 7.0
  x <- 4.0
  for (transition in transitions) {
    if (transition %in% names(commands)) {
      # Get the command and properly escape it
      command <- commands[transition]
      
      # Handle escaped quotes properly in XML
      command <- gsub('"', '&quot;', command)
      delay_attr <- sprintf('delay="%s"', command)
    } else {
      # Default delay attribute if no command
      delay_attr <- 'delay="1.0"'
    }
    
    transitions_xml <- paste0(transitions_xml, 
                              sprintf('\n      <transition %s delay-x="0.0" delay-y="0.0" label-x="0.0" label-y="0.0" name="%s" nservers-x="0.5" rotation="0.0" type="EXP" x="%s" y="%s"/>', 
                                      delay_attr, transition, x, y))
    x <- x + 3.0
    if (x > 15) {
      x <- 4.0
      y <- y + 3.0
    }
  }
  
  # Edges/Arcs
  edges_header <- '\n    </nodes>\n    <edges>'
  edges_xml <- ""
  
  for (i in 1:nrow(arc_df)) {
    row <- arc_df[i, ]
    mult_attr <- ''
    if (!is.na(row$multiplicity) && row$multiplicity > 1) {
      mult_attr <- sprintf(' mult="%d"', row$multiplicity)
    }
    
    if (row$direction == "INPUT") {
      edges_xml <- paste0(edges_xml, 
                          sprintf('\n      <arc head="%s" kind="INPUT"%s tail="%s"/>', 
                                  row$transition, mult_attr, row$place))
    } else if (row$direction == "OUTPUT") {
      edges_xml <- paste0(edges_xml, 
                          sprintf('\n      <arc head="%s" kind="OUTPUT"%s tail="%s"/>', 
                                  row$place, mult_attr, row$transition))
    }
  }
  
  # Closing tags and standard measures
  xml_footer <- '\n    </edges>
  </gspn>
  <measures gspn-name="PetriNet" name="Measures" simplified-UI="false">
    <assignments/>
    <greatspn/>
    <formulas>
      <formula comment="Basic statistics of the toolchain execution." language="STAT"/>
      <formula comment="All the basic Petri net measures" language="ALL"/>
    </formulas>
  </measures>
</project>'
  
  # Assemble the complete XML
  xml_content <- paste0(xml_header, places_xml, transitions_xml, edges_header, edges_xml, xml_footer)
  
  # Write to file
  tryCatch({
    # Ensure the output directory exists
    dir.create(dirname(pnpro_out), showWarnings = FALSE, recursive = TRUE)
    
    # Write the XML content to file
    writeLines(xml_content, con = pnpro_out)
    message("Successfully wrote PNPRO file to: ", pnpro_out)
  }, error = function(e) {
    stop("Error writing PNPRO file: ", e$message)
  })
  
  return(pnpro_out)
}