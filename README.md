# Detecting-Face-and-Torso-
With this Code , I found out the possible Head and Torso from Sports Images Separately and connected them all by matching their orientation , so if Person x 's Face is Detected it will not map into Person y's Torso , it will map into Torso of x only irrespective of Torso Not Getting Detected , so with 100 % Accuracy , i detected Face and Torso

Here I have faced 4 issues / cases :
a. Face and Body Detected
	b. Face Detected but Body not Detected
	c. Face not Detected but Body Detected
	d. Face and Body not Detected

a. Face and Body Detected :  For this case we have no issue we will move to next step
b.Face Detected but Body not Detected: For this we consulted a paper where it was mentioned height of body upto thighs is 7 times the height of face , so we expanded the face downwards towards opposite boundary upto 7*height of face length
c.Face not Detected but Body Detected: Same as above we will move the body upwards towards opposite boundary until my total height is total height / 7
d.Face and Body not Detected: Here we use Boundary Growing Method of detecting Undetected face and then expanding towards down boundary to get the undetected body
