This project was produced for PHS866 as part of the Clinical Scientific Computing course at Liverpool University. 

This document is a readme to aid in using the project the web application is currently live at https://phys866-model-sam-mccaughran-513d18b2c906.herokuapp.com/

The files for this project are also available at https://github.com/Sam-McCaughran/PHYS866_Part1

<h2>Structure</h2>
<h3>File List</h3>
The file includes a number of .png files, csv files, jupyter notebooks and a python file - these are:

<ul>
    <li>ModelParameters.png is the parameters used for my 'noisy' model.</li>
    <li>Papers1.png is the full paper list from the meta analysis.</li>
    <li>Papers2.png is the paper list used for collating values.</li>
    <li>Presentation.png was the background used during the presentation.</li>
    <li>1_PaperAnalysis.ipynb is the script used to analyse the data from the meta analysis.</li>
    <li>streamlit_app.py is the source code for the web application.</li>
    <li>3_Data_Analysis.ipynb is the analysis file for experimental results.</li>
    <li>Paper_Analysis_Data.csv is the output of the 1_PaperAnalysis.ipynb.</li>
    <li>AnalysisResults.csv is the output of 3_Data_Analysis.ipynb.</li>
</ul>

<h3>Sub Folders</h3>
<ul>
    <li>Figures holds all of the output figures from 3_Data_Analysis.ipynb.</li>
    <li>Experimental data is where any data exported from the web app needs to be place to be automatically picked up by 3_Data_Analysis.ipynb.</li>
</ul>
<h2>Run Web App and Scripts Locally</h2>

To run the model locally please follow the steps below:
<ol>
    <li>Ensure all dependencies are installed.</li>
    <li>Run streamlit_app.py using the command streamlit run streamlit_app.py </li>
    <li>Navigate to the address provided by the terminal and create experimental data.</li>
    <li>Export the data using the 'Export Data Button'.</li>
    <li>Place the exported csv in the 'Experimental_Data' folder - note the file must contain 'cell_survival_experiments' with anything else appended - delete the existing data if you do not want it to be included in your analysis.</li>
    <li>Run 3_Data_Analysis.ipynb -  At present the plots in this file only works for liver values as I did not have time to make general functions however this would be relatively easy to achieve.</li>
</ol>
If you have any questions please don't hestiate to contact me at sam.mccaughran@gmail.com
