Assignment #3: Ray tracing



FULL NAME: Ching-Chih Chen




MANDATORY FEATURES

------------------



<Under "Status" please indicate whether it has been implemented and is

functioning correctly.  If not, please explain the current status.>



Feature:                                 Status: finish? (yes/no)

-------------------------------------    -------------------------

1) Ray tracing triangles                  !!!yes!!!


2) Ray tracing sphere                     !!!yes!!!


3) Triangle Phong Shading                 !!!yes!!!


4) Sphere Phong Shading                   !!!yes!!!


5) Shadows rays                           !!!yes!!!


6) Still images                           !!!yes!!!
   

7) Extra Credit (up to 10 points)
   
    !!! Recursive reflection: ray will keep reflecting until ks1*ks2*... < 0.01 or no 
	intersection, and last reflection light will be ambient light. Pixel color will be
	(1-ks1) * localPhongColor + ks1 * ((1-ks2) * localPhongColor2 + ks2 * (...)).!!!

    !!! antialiasing: For every pixel, I use average of 4 rays to render, like double image
	size and averaging every 4 pixels !!!
   
    !!! Soft shadows: I turn the point light into square area light with 
	different x,y and same z values, and each light intensity will
	become ( 1/lightarea ) !!!
