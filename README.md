# Upload-Application

For an in depth description of the Pseudotime code, please see 'Overview of Pseudotime Analysis Application.pdf'. The necessary packages can be seen in sessionInfo(). A copy of the thesis associated with this project is available here.

'tcell_libs_raw_dash.R' contains the functions to produce the objects to be loaded by the application ('app.R'). These objects are not pre-processed and are created by the user to be saved in the specified folder (working directory). See ' Upload Application Steps.pdf' for amore in-depth look at the steps to run the Upload Application. 

The working directory should be changed in the code to load/save files (setwd()). It can be done on line 4 of app.R
![image](https://user-images.githubusercontent.com/88155528/138607672-3ba9c5a0-0898-44b1-8354-cb6eace55af4.png)

Please ensure the necessary packages and versions are installed. If issues, occur with PhateR see 'Upload Application Steps.pdf' on the steps to take.

![image](https://user-images.githubusercontent.com/88155528/138607713-838fe26b-57c1-4812-a598-3b0a65f769e3.png)
![image](https://user-images.githubusercontent.com/88155528/138607752-51f68258-b870-4827-91b7-f6df4ad42f58.png)


The Seurat files could not be uploaded as they are too big and could not be compressed. Use these to produce Seurat objects if you wish to test the application. If you wish to test with resolution 0.15 set both the highest and lowest slider to the resolution. It should take around 1.5 hours but may vary depending on the device
WT1.data <- Read10X(data.dir = "/data/2623287c/Project1/newtcell/t_cell_raw_data/WT1")
WT1 <- CreateSeuratObject(counts = WT1.data, project = "WT1", min.cells = 3)
