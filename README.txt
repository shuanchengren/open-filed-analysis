#Title of code: Code for open field analysis for the article “A midbrain GABAergic circuit constrains wakefulness in stressful challenging conditions”

*Features: MATLAB code for tracking animals’ locomotor activity in the open filed test. This code permits tracking the trajectory of animals in the open field and analyzing the travel speed, distance traveled. 
*Format(s): .mat.
*Size(s): 72kB.
Code was run on MATLAB (2016b) on Windows 10. “MY_MICE_TRACKING_56_GENERALSIZE” is code for generating “my_mice_tracking_56_generalSize.fig”. MY_MICE_TRACKING_56_GENERALSIZE, by itself, creates a new MY_MICE_TRACKING_56_GENERALSIZE or raises the existing. Detailed instructions can be found in the scripts. 
Load: 
**Tutorials**
1. Press the “Load AVI” button to load an AVI format video file;
2. Randomly select one video picture by dragging the scroll bar under the Main Window. Press the purple “Set Mask” button, the video picture in the Main Window will the previewed in the right small Preview Window;
3. Randomly select another different video picture by dragging the scroll bar under the Main Window. Please make sure that animal is not in the same position in these two pictures;
4. Press the “Set Mask” button again, the moving animal will be recognized and the disappeared in the Preview Window. Then, moving the scroll bar back to the start. 
5. Press the cyan “Mice tra…” button, colored picture in the Main Window will be transformed into gray scale picture and the animal will be identified as white. In addition, a “Select a rectangle to hold the mouse” window will be presented;
6. Press the “确定” button and select a rectangle area to cover the area of mouse in the Main Window. The animal will be zoomed in the Preview Window. A “Setup polygon to hold the mouse” window will be presented;
7. Press the “确定” button and draw the outline the animal in the Preview Window by click the left mouse button. Double click left mouse button, the outline of the animal in the Preview Window will disappear and the tracking process with start. The position of animal tracking by every frame and will be showed in the Main Window in every 500 frames. 
8. Please note there is a “Error Check” button under the “Mice tra…” button. You can enter a number to set the threshold for checking the tracking error. The number indicates the pixel changes in two adjacent video frames. The tracking error will be presented under the main window. You can check the error by locating to the specific frame.
9. After finish tracking the location of animal, the “AVI Parameter Setup Panel” can be used to calculate the traveled pixels of animal. Entering the “Start” and “End” frame and then pressing the “Calculate Path” button to get the length of path. 
10. Left to the “AVI Parameter Setup Panel”, “Rect Plot” button is used to generate the traveled path plot and the “TimeM” button is used to generate the heatmap of time of animal spend in each location of the experimental area. 


