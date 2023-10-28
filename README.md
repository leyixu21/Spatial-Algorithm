# Spatial-Algorithm
The course project Spatial Algorithm is a group project with three group members. This project utilized crowdsourced images to detect cultural ecosystem services in the canton of Graubünden, Switzerland with Python, and it is a replicate of the paper by Gianfranco Gliozzo [[1]](#1). The algorithms in this project were implemented manually without relying on existing Python packages.

## First Algorithm: Rasterize Polygon
The goal of this algorithm, rasterize polygon, was, as the name indicates, to rasterize the polygon that was created from the border of the Canton of Graubünden. The algorithm requires three inputs: 
1) A list of centroid coordinates of all raster cells within the bounding box of the border polygon
2) The resolution of the raster, either 1x1km or 2x2km
3) The polygon of the border of the Canton of Graubünden that is defined as a Polygon class

The output is a list only containing two separate values: “1” or “None”. A “1” indicates that the raster cell is inside the polygon. Simultaneously, “None” indicates that the raster cell is outside the polygon and thus outside of the Canton of Graubünden. The length of the output list is identical to the list of centroid coordinates since every raster cell gets one value.

## Second Algorithm: Point Index
This is the longest algorithm in our code and consists of almost 200 lines. The algorithm builds a spatial index for the raster cells within the bounding box of the Canton of Graubünden with a specific resolution. It counts the number of points (contributors in this case) per raster. The input of this algorithm includes:
1) The images as a Point Class
2) Their corresponding user ID, to get the contributors, not the images
3) The bounding box of the Graubünden polygon
4) The resolution of raster cells

Based on the spatial index of input (1), two lists are created for each raster: one list stores the user ID, and the other one stores the number of contributors. If the user ID is already in the raster cell, it won’t be included again. The second list will be used to calculate Moran’s I.

## Third Algorithm: Moran’s I Index
This algorithm calculates the global Moran’s I, the local Moran’s I, and the expected value of Moran’s I. The input is the nCols*nRows list that was computed in the second algorithm and stores the value for each raster cell. In this case, it should be the reshaped list storing the number of contributors per raster cell. To calculate Moran’s I, the queen weight matrix was applied. Three functions were globalMoranI, localMoranI, and expectedValue.

## Examples
- Concentration of image contributors per cell, 2x2 km grid resolution
<img src="https://github.com/leyixu21/Spatial-Algorithm/assets/96665869/76078375-4133-4666-a746-d7236986652c" height="300">

- Local Moran's I, 2x2 km grid resolution
<img src="https://github.com/leyixu21/Spatial-Algorithm/assets/96665869/96135862-ca5c-406c-8dcc-35d16257968f" height="300">

## References
<a id="1">[1]</a> 
G. Gliozzo, N. Pettorelli, and M. Muki Haklay, “Using crowdsourced imagery to detect cultural 
ecosystem services: A case study in South Wales, UK,” Ecology and Society, vol. 21, no. 3, 2016, 
doi: 10.5751/ES-08436-210306.
