
First we need to find the string point of the trajectory, we use the 
car's current s and d value as starting point if the previous_path has 
less than 2 elements, otherwise we use the last two points in the 
previous_path in order to make the path smooth.

Then we use current s position plus 30, 60 and 90 as s position on
the path, d will depend on whether we will keep lane or change 
lane left or right. 

Now we have 3 points in Frenet coordinate system, and we convert
them to map coordinate system using getXY, then we convert the points
from map coordinate system to ego vehicle local coordinate 
system. 

We feed the three points to spline and get the result.

We add the points in the previous_path to the list of points
that will be sent to the simulator.

We giving x values in the ego vehicle coordinate system
 based on the reference speed to the fitted 
spline curve and get y values.

We convert these points to the map coordinate system and add 
them to the list of points that will be sent.