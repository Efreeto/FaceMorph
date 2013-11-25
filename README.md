FaceMorph
=========

Takes two 3D head models and make them morph to each other. (June 2010)

------------------------------------------------------------------

Usage:  
Left click and drag - view change  
Right click - open up menu to select the second head (morphs to the first one)  
'p' key - pause the demo  
Space bar - show/hide control points  

It's a simple code sample I uploaded to show companies. I utilized a few libraries for math, visualization, and most importantly, the RBF interpolation algorithm. On the launch, a head starts morphing into a different head and when it does, the head shape comes back to its original. Facial expressions can also be achieved by changing the control point locations.  
If you want to make it work for your own head models, you must provide their corresponding vertex index pairs or provide new locations.
