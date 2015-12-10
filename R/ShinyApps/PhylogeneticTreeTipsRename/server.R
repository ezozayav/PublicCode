#author: Mark B Schultz
#email: dr.mark.schultz@gmail.com
#date: 04 Sep 2015

library(rncl)
library(ggtree)
library(ggplot2)
library(ape)

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 9MB.
options(shiny.maxRequestSize = 20*1024^2)

shinyServer(function(input, output) {
    
    # reactive to bring up the input tree file
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    # the tree object will be stored in "treeInput()"
    treeInput <- reactive({
        inFile = input$file1
        if (is.null(inFile))
            return(NULL)
        if(input$tree_informat=="newick"){
            tree = read_newick_phylo(inFile$datapath)
            return(tree)
            }
        if(input$tree_informat=="nexus"){
            tree = read_nexus_phylo(inFile$datapath)
            return(tree)
            }
        if(input$tree_informat=="beast"){
            beast_tree = read.beast(inFile$datapath)
            return(beast_tree)
            }
        })
    
    # reactive to bring up the input metadata file
    # the metadata will be stored in "metadataInput()"
    metadataInput <- reactive({
        inFile2 = input$file2
        if (is.null(inFile2))
            return(NULL)
        table = read.csv(inFile2$datapath, header=input$header, sep=input$sep, 
                         quote=input$quote, as.is = TRUE, check.names = FALSE)
        return(table)
    })
    
    # interactive to bring up the drop down to choose original strain names 
    # column after loading input metadataInput()
    output$choose_original_label_column = renderUI({
        if (is.null(metadataInput()))
            return(NULL)
        colnames = names(metadataInput())
        selectInput("old_names", "3) Metadata column containing original tip labels:", as.list(colnames))
    })
    
    # interactive to choose column headers of the input metadata to paste 
    # together to make new tip names
    output$choose_new_label_columns = renderUI({
        if (is.null(metadataInput()))
            return(NULL)
        colnames = names(metadataInput())
        checkboxGroupInput("new_names", "New names (select variables to paste together):", choices=colnames, selected="")
    })

    output$tree_in = renderPlot({
        if (is.null(treeInput()))
            return(NULL)
        if(input$tree_informat == "newick" | input$tree_informat == "nexus"){
            plot(treeInput())}
        if(input$tree_informat == "beast"){
            ggtree(treeInput(), showDistance = TRUE)  + 
            geom_tiplab(size=3, align=FALSE) +
            theme_tree2()
            }
    })

    output$metadata = renderTable({
        if (is.null(metadataInput()))
            return(NULL)
        metadataInput()
    })
    
    # rename tips function newick and nexus (and beast tree@phylo element)  
    rename_tips = function(tre, metadata_table)
    {
        for(i in 1:length(tre$tip.label))
        {
            new_tip = paste0(as.character(metadata_table[which(metadata_table[,input$old_names]==tre$tip.label[i]),c(input$new_names)]), sep="", collapse="_")
            tre$tip.label[i] = as.character(new_tip)
        }
        return(tre)
    }
    
    # rename translation column in beast tree@translation matrix
    rename_trans_matrix_col3 = function(trans_matrix, metadata_table)
    {
        trans_rename_col3 = c()
        for(i in 1:length(trans_matrix[,2]))
        {
            trans_rename_each = paste0(as.character(metadata_table[which(metadata_table[,input$old_names]==trans_matrix[,2][i]),c(input$new_names)]), sep="", collapse="_")
            trans_rename_col3 = c(trans_rename_col3, as.character(trans_rename_each))
        }
        return(trans_rename_col3)
    }

    output$tree_out = renderPlot({
        if (is.null(treeInput()))
            return(NULL)

        #plot the tree, newick or nexus
        if(input$tree_informat == "newick" | input$tree_informat == "nexus"){
            tree_renamed = rename_tips(treeInput(), metadataInput())
            plot(tree_renamed)
        }
        if(input$tree_informat == "beast"){
            # store the tree locally
            tree_local = treeInput()
            # rename the tips in the local tree
            tree_renamed = rename_tips(tree_local@phylo, metadataInput())
            # add node labels, e.g. Node1, Node2..NodeX to the renamed tree
            tree_renamed = makeNodeLabel(tree_renamed)
            # reverse the node.labels
            tree_renamed$node.label = rev(tree_renamed$node.label)
            nodelabels_tree_renamed = tree_renamed$node.label
            # calculate what the node number would be for Node1
            node_1_node_number = length(tree_renamed$tip.label)+1
            # calculate what the node number would be for NodeX
            node_n_node_number = nrow(tree_renamed$edge)+1
            # create a translation table between Node1..NodeX and node numbers 
            node_numbers = c(node_1_node_number:node_n_node_number)
            node_trans =  cbind(label=nodelabels_tree_renamed, number=node_numbers)
            # store the tree name to node number translation matrix locally
            tree_trans = treeInput()@translation
            # store the beast tree node metadata locally
            data = treeInput()@stats
            # get the renamed tree tips for cbinding to the local translation matrix
            trans_matrix_col3 = rename_trans_matrix_col3(tree_trans, metadataInput())
            # bind the renamed tips
            tree_trans = cbind(tree_trans, renamed=trans_matrix_col3)
            # add a label column to the start of the locally stored beast node metadata
            data = cbind(label=rownames(data), data)
            # remove the factor formatting from the newly added column
            data[,"label"] = as.character(data[,"label"])
            # flatten the lists in the data cols as many of the cols are lists of lists
            for(i in 1:ncol(data)){
                for(j in 1:nrow(data)){
                    data[j,i] = toString(unlist(data[j,i], use.names=F))
                    data[j,i] = gsub("NA", "", data[j,i])
                }
                data[,i] = unlist(data[,i], use.names=F)
            }
            # populate the data[,"label"] column with the node.labels from node_trans
            for(i in 1:nrow(data)){
                for(j in 1:nrow(node_trans)){
                    if(data[i,"label"] == node_trans[j,"number"]){
                        data[i,"label"] = node_trans[j,"label"]
                    }
                }
            }
            # populate the data[,"label] column with the translated tree tips after doing tree_renamed
            for(i in 1:nrow(data)){
                if(as.numeric(rownames(data[i,]))<=nrow(tree_trans)){
                    a = rownames(data[i,])
                    data[i,"label"] = toString(tree_trans[which(tree_trans[,1]==a),3])
                    tree_trans[which(treeInput()@translation[,1]==a), 3]
                }
            }

            # remove the node column at the end of the matrix
            data[,"node"] = NULL
 
            # ladderize the tree, ascending
            tree_renamed = ladderize(tree_renamed, right = FALSE)
            # store the tree in a tempfile() connection
            tmpFile = tempfile()
            write.tree(tree_renamed, file=tmpFile)
            outfile <- tempfile()
            write.jplace(tmpFile, data, outfile)
            # read in the tree from the tempfile 
            jp <- read.jplace(outfile)
            # plot the tree stored in jp
            plot(jp@phylo, cex=0.8)
        }
        
        # specify the file extension for download
        if(input$tree_informat == "newick" | input$tree_informat == "nexus"){
            extension = "_tips_renamed_nwk.tre"}
        if(input$tree_informat == "beast"){
            extension = "_tips_renamed_JPLACE.tre"}
        
        # create the download button and data
        output$downloadData <- downloadHandler(
            filename = function() {
                paste(input$file1, extension, sep='')
            },
            content = function(file) {
                if(input$tree_informat == "newick" | input$tree_informat == "nexus"){
                    write.tree(tree_renamed, file)}
                if(input$tree_informat == "beast"){
                    write.jplace(tmpFile, data, file)
                }
            }
        ) 
    # end brace for output$tree_out renderPlot tab
     })
    
    # The following text will be plotted at the top of the main panel if there 
    # are metadata columns selected
    output$text1 = renderText({
        if(length(input$new_names)>0){
        print(paste("Pasting", length(input$new_names), "columns together to make new tip names", sep= " "))
        }
    })
})
