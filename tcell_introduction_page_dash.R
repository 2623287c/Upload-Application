HTML("
     
<div>
  <h1>Welcome</h1>
  <p>
  Here, you can explore the single-cell RNA sequencing (scRNA-Seq) results from Seurat objects you upload. The application typically can take between 1.5-6 hours per resolution, depending on the number of cells. Please change the working directory (line 4 of app.R) to your working directory to ensure you can load/save files. Please see GitHub (https://github.com/2623287c/Upload-Application) for further details
  <p>
  
  <p>
  The analysis was carried out using the Seurat pipeline for differential expression. It uses the integrated workflow so requires more than one group to integrate. 
  </p>
  
  </p>
  
  
  <div>
    <b>Clustering</B> - Allows cluster visualization and exploration of top cluster markers. Can change resolution and view the results. Can also change cluster names (note: cluster names will not be updated in the Pseudotime section).
    </div>
    
    
    <div>
     <b>Differential expression (DE)</b> - Comparison of gene expression.  
    
  
    </div>
    
      
    <div>
     <b>Pseudotime </b> - Four different tools (Slingshot, tradeSeq (downstream of Slingshot), Monocle 2 and Monocle 3) allowing pseudotime analysis of the data. Plots of the trajectories and heatmaps are available for each tool. The method to create the heatmaps varies for each tool and are not directly comparable, all except tradeSeq (due to the nature of the tool) are ordered by qvalue with the top 40 genes shown in the heatmap. The tradeSeq, Monocle 2 and Monocle 3 heatmaps can take between 5 and 10 minutes to load.
  
    </div>
  
</div>")