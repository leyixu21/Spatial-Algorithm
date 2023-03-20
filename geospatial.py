# Classes and methods for geospatial algorithms based on points with x, y coordinates

# for GEO877 Spatial Algorithms
# - requires Python 3.6 or later (or replace f-string with older print format)
# - Point distance measures return results in metres (or assume metres)

# release history:
# 2021-03-27, v1.0: updated version of points.py (v1.0) for solution 1 (comparing distance measures)
# author:  Sharon Richardson
# version: 1.0

# 2022-03-25, v2.1: Modified Sharon's original polygon code slightly
# author:  Ross Purves


# 2022-04-02, v2.2: Added Segment class and determinants to Point class
# author:  Ross Purves
# version: 2.2
# --------------------------
# 2022-06-04, v3.0:
# 1. Modified PointIndex class to count the number of contributors per raster
# 2. Added borderRaster class to rasterize the border polygon
# 3. Added MoranI class to calculate the global Moran's I, local Moran's I, and expected value
# --------------------------

from re import X
from tkinter import E
from numpy import sqrt, radians, arcsin, sin, cos, square

# class and methods for a geometric point
# =======================================
from numpy import sqrt
from pandas import NA


####### Point #######

class Point():
    # initialise
    def __init__(self, x=None, y=None):
        self.x = x
        self.y = y
    
    # representation
    def __repr__(self):
        return f'Point(x={self.x:.2f}, y={self.y:.2f})'

        # Test for equality between Points
    def __eq__(self, other): 
        if not isinstance(other, Point):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return self.x == other.x and self.y == other.y

    # We need this method so that the class will behave sensibly in sets and dictionaries
    def __hash__(self):
        return hash((self.x, self.y))
    
    # calculate Euclidean distance between two points
    def distEuclidean(self, other):
        return sqrt((self.x-other.x)**2 + (self.y-other.y)**2)
    
    # calculate Manhattan distance between two points
    def distManhattan(self, other):
        return abs(self.x-other.x) + abs(self.y-other.y)

    # Haversine distance between two points on a sphere - requires lat/lng converted to radians
    def distHaversine(self, other):
        r = 6371000  # Earth's radius in metres (will return result in metres)
        phi1 = radians(self.y) # latitudes
        phi2 = radians(other.y)
        lam1 = radians(self.x) # longitudes
        lam2 = radians(other.x)

        d = 2 * r * arcsin(sqrt(sin((phi2 - phi1)/2)**2 + 
                                      cos(phi1) * cos(phi2) * sin((lam2 - lam1)/2)**2))
        return d   

    
    # Calculate determinant with respect to three points. Note the order matters here - we use it to work out left/ right in the next method
    def __det(self, p1, p2):
        det = (self.x-p1.x)*(p2.y-p1.y)-(p2.x-p1.x)*(self.y-p1.y)       
        return det

    def leftRight(self, p1, p2):
        # based on GIS Algorithms, Ch2, p11-12, by Ningchuan Xiao, publ. 2016
        # -ve: this point is on the left side of a line connecting p1 and p2
        #   0: this point is collinear
        # +ve: this point is on the right side of the line
        side = int(self.__det(p1, p2))
        if side != 0:
            side = side/abs(side)  # will return 0 if collinear, -1 for left, 1 for right
        return side

####### Segment #######
class Segment():    
    # initialise
    def __init__(self, p0, p1):
        self.start = p0
        self.end = p1
        self.length = p0.distEuclidean(p1)
    
    # representation
    def __repr__(self):
        return f'Segment with start {self.start} and end {self.end}.' 

    # Test for equality between Segments - we treat segments going in opposite directions as equal here
    def __eq__(self, other): 
        if (self.start == other.start or self.start == other.end) and (self.end == other.end or self.end == other.start):
            return True
        else:
            return False
            # We need this method so that the class will behave sensibly in sets and dictionaries
    
    def __hash__(self):
        return hash((start.start.x, self.start.y,self.end.x, self.end.y))    
        
    # determine if intersects with another segment (using Point method leftRight)
    # - should we incorporate testing for identical segments and non-zero lengths?
    def intersects(self, other):        
        # create bounding boxes for each segment, using Bbox class
        self_bbox = Bbox(self)
        other_bbox = Bbox(other)

        # test if bounding boxes overlap (using Bbox method testOverlap)
        bbox_overlap = self_bbox.intersects(other_bbox)
               
        # if the two bboxes do not overlap, there can be no intersection        
        if bbox_overlap == False:
            return False
        
        else: # bboxes overlap, test if lines intersect or are collinear
            
            # find side for each point of each line in relation to the points of the other line
            apq = self.start.leftRight(other.start, other.end)
            bpq = self.end.leftRight(other.start, other.end)
            pab = other.start.leftRight(self.start, self.end)
            qab = other.end.leftRight(self.start, self.end)           
            
            # leftRight sums to 0 for both segments if they intersect (+1/-1) or are collinear (0/0)
            if (apq + bpq == 0 and pab + qab == 0): 
                return True
            else:
                return False      
    
####### Bounding Box #######    

# define bounding box for 2 or more points (initialised as PointGroup, Polygon or Segment)
class Bbox():
    
    # initialise
    def __init__(self, data):
    # using built-in `isinstance` to test what class has been used to initialise the object   

        # for Segment objects 
        if isinstance(data, Segment) == True:
            x = [data.start.x, data.end.x]
            y = [data.start.y, data.end.y]
        
        # for PointGroup objects (including Polygon)
        else:      
            x = [i.x for i in data]   # extract all x coords as a list
            y = [i.y for i in data]   # extract all y coords as a list

        # determine corners, calculate centre and area
        self.ll = Point(min(x), min(y))    # lower-left corner (min x, min y)
        self.ur = Point(max(x), max(y))    # upper-right corner (max x, max y)
        self.ctr = Point((max(x)-min(x))/2, (max(y)-min(y))/2)   # centre of box
        self.area = (abs(max(x)-min(x)))*abs((max(y)-min(y)))    # area of box
           
    # representation
    def __repr__(self):
        return f'Bounding box with lower-left {self.ll} and upper-right {self.ur}' 
    
    def __eq__(self, other): 
        if (self.ll == other.ll and self.ur == other.ur):
            return True
        else:
            return False
            # We need this method so that the class will behave sensibly in sets and dictionaries
    
    def __hash__(self):
        return hash(self.ll, other.ur)  
        
    # test for overlap between two bounding boxes
    def intersects(self, other):       
        # test if bboxes overlap (touching is not enough to be compatible with the approach to segments)
        if (self.ur.x > other.ll.x and other.ur.x > self.ll.x and
            self.ur.y > other.ll.y and other.ur.y > self.ll.y):
            return True
        else:
            return False
        
    # Check if bounding box is contained by this one     
    def contains(self, other):
        if (self.ur.x >= other.ur.x and self.ll.x <= other.ll.x and
            self.ur.y >= other.ur.y and self.ll.y <= other.ll.y):
            return True
        else:
            return False
        
    def intersectsRegion(self, other):
        if self == other:
            return self
        elif self.contains(other):
            return other
        elif self.intersects(other):
            llx = max(self.ll.x, other.ll.x)
            lly = max(self.ll.y, other.ll.y)
            urx = min(self.ur.x, other.ur.x)
            ury = min(self.ur.y, other.ur.y)
            return Bbox([Point(llx, lly), Point(urx, ury)])
        else:
            return None
                        
    # contains includes points on boundaries, otherwise we have problems when points are used to define
    # box but not in it
    def containsPoint(self, p):
        if (self.ur.x >= p.x and p.x >= self.ll.x and
            self.ur.y >= p.y and p.y >= self.ll.y):
            return True
        return False
        

####### Point Group #######    

# FOR GROUPS OF POINTS

# data provided should be as an array of points with x, y coordinates.

# class for a group of Points, assumes initial data is unsorted, spatially
class PointGroup(): 
    # initialise
    def __init__(self, data=None, xcol=None, ycol=None):
        self.points = []
        self.size = len(data)
        for d in data:
            self.points.append(Point(d[xcol], d[ycol]))
    
    # representation
    def __repr__(self):
        return f'PointGroup containing {self.size} points' 
 
    # create index of points in group for referencing
    def __getitem__(self, key):
        return self.points[key]
    
####### Polygon #######
# Polygon class for polygons, assumes initial data is in a spatially sorted order
class Polygon(PointGroup):  
    # initialise
    def __init__(self, data=None, xcol=None, ycol=None):
        self.points = []
        self.size = len(data)
        for d in data:
            self.points.append(Point(d[xcol], d[ycol]))
        self.bbox = Bbox(self)
        
    # representation
    def __repr__(self):
        return f'Polygon PointGroup containing {self.size} points' 
  
    # test if polygon is closed: first and last point should be identical
    def isClosed(self):
        start = self.points[0]
        end = self.points[-1]
        return start == end

    def removeDuplicates(self):
        oldn = len(self.points)
        self.points = list(dict.fromkeys(self.points)) # Get rid of the duplicates
        self.points.append(self.points[0]) # Our polygon must have one duplicate - we put it back now
        n = len(self.points)
        print(f'The old polygon had {oldn} points, now we only have {n}.')
        
        # find area and centre of the polygon
    # - based on GIS Algorithms, Ch.2 p9-10, by Ningchuan Xiao, publ. 2016 
    
    def __signedArea(self):  # used for both area and centre calculations - this is a private method (only used within the class)
        a = 0
        xmean = 0
        ymean = 0
        for i in range(0, self.size-1):
            ai = self[i].x * self[i+1].y - self[i+1].x * self[i].y
            a += ai
            xmean += (self[i+1].x + self[i].x) * ai
            ymean += (self[i+1].y + self[i].y) * ai

        a = a/2.0   # signed area of polygon (can be a negative)
    
        return a, xmean, ymean
    
    def area(self):
        a, xmean, ymean = self.__signedArea()
        area = abs(a)   # absolute area of polygon
        
        return area

    def centre(self):
        a, xmean, ymean = self.__signedArea() # note we use the signed area here
        centre = Point(xmean/(6*a), ymean/(6*a)) # centre of polygon 
        return centre

    def containsPoint(self, p):
        if (self.bbox.containsPoint(p) == False):
            return False
        
        # Solution, as discussed in lecture, added here
        ray = Segment(p, Point(self.bbox.ur.x+1, p.y))
        count = 0
        
        for i in range(0, self.size-1):
            start = self[i]
            end = self[i+1]
            if (start.y != end.y): 
                if ((p.x < start.x) and (p.y == start.y)): 
                    count = count + 1
                else:
                    s = Segment(start, end)
                    if s.intersects(ray):
                        count = count + 1
              
        if (count%2 == 0):
            return False           

        return True

####### Point Index #######
class PointIndex():
    # initialise the index
    def __init__(self, data, user, box, res):                        
        # Work out number of rows and colums
        self.res = res
        w = box.ur.x - box.ll.x
        h = box.ur.y - box.ll.y
        self.nCols = int(w/ self.res) + 1 # +1 is to make sure the grid covers all of our points
        self.nRows = int(h/ self.res) + 1
        self.points_num = [] # Store the number of contributors in each raster      
        self.points = [] # Store the user id in each raster
        self.raster_coords = [] # store the centroid coordinate of each raster
        
        # Bounding box should be redefined to be same as grid now; we will use it later
        ur = Point(box.ll.x + (self.nCols * self.res), box.ll.y + (self.nRows * self.res))
        ll = box.ll
        self.bBox = Bbox([ll, ur])
        self.maxIndex = (self.nCols * self.nRows) -1 # Store the maximum possible index value, to use as a check
        self.bigArray = [] # This array is used to carry out brute force calcuations; it should be removed
        # Fill an empty index
        i = 0
        while i <= self.maxIndex:
            self.points_num.append(0) # Store the number of contributors
            self.points.append([]) # Store a list of contributor id
            i += 1
            
        self.addPoints(data, user)
        self.raster_coords = self.__getCoords()


    # Pretty print of the data structure - note how we make sure things are displayed in the right order    
    def __repr__(self):
        s = f'Resolution: {self.res}\n'
        s = s + f'nCols: {self.nCols}\n'
        s = s + f'nRows: {self.nRows}\n'
        count = 0
        for i in range(1, self.nRows+1):
            for c in range(0,self.nCols):                
                r = self.nRows - i
                #print(f'{r} {c}')
                index = (r * self.nCols) + c
                s = s + f'{self.points_num[index]}\t'
                count = count + self.points_num[index]            
            s = s + '\n'
        s = s + f'Points: {count}'
        return s
    
    # Method to add an array of point data - can add more data later, IF it is in our bounding box 
    def addPoints(self, data, user):
        for p, u in zip(data, user):
            self.addPoint(p, u)
            self.bigArray.append(p) # Single bucket with all points

    # Put an individual point in the correct bucket, and increment its count        
    def addPoint(self, p, u):
        i = self.pointIndex(p)
        if self.points_num[i] == 0:
            self.points_num[i] = 1
            self.points[i] = [u]
        else:
            # If the user id is not in raster i, count it
            if u not in self.points[i]:
                count = self.points_num[i]
                ps = self.points[i]
                ps.append(u)
                count +=1
                self.points_num[i] = count
                self.points[i] = ps
        
    # Get back the index of a point using its coordinates
    def pointIndex(self, p):
        j = int((p.y - self.bBox.ll.y)/ self.res)
        i = int((p.x - self.bBox.ll.x)/ self.res)
        return (j * self.nCols) + i  

    # Check if a point is stored in the index
    def contains(self, p):
        i = self.pointIndex(p)
        ps = self.points[i][1]
        return p in self.points
    
    # Get back all the points in a region
    def regionQuery(self, region):        
        # We first return the query region that intersects with our bounding box
        query = self.bBox.intersectsRegion(region)

        # If the query is wholly outside our index we return an empty list
        if query is None:
            return [], 0
        
        # Define the bounds of our search for points
        xStop = query.ur.x + self.res
        yStop = query.ur.y + self.res
        ps = []

        x = query.ll.x
        while (x <= xStop):
            y = query.ll.y
            while (y <= yStop):
                i = self.pointIndex(Point(x,y))
                if (i < self.maxIndex and self.points[i][0] > 0): # If bucket lies in index, and it contains points then
                    ps = [*ps, *self.points[i][1]]                # merge the points from this bucket
                y = y + self.res
            x = x + self.res
        
        #Now iterate through the points in the buckets, and only keep those in our query region
        final = []
        for p in ps:
            if (region.containsPoint(p)):
                final.append(p)
                
        return final, len(final)
                
    def nearestPoint(self, p):
        index = self.pointIndex(p)
        ps = []
        if (index >= 0 and index <= self.maxIndex):
            ps = self.points[index][1] # Get the points in this cell
        # If there are no points in the cell we look for more
        # Increase the range we search by the cell resolution
        step = self.res
        while len(ps) == 0:
            ll = Point(p.x - step, p.y - step)
            ur = Point(p.x + step, p.y + step)
            box = Bbox([ll,ur])
            ps, count = self.regionQuery(box)
            #print(f'Loop \n{p} \n{box}  \n{ps}')
            step = step + self.res

        oldCount = len(ps)
        #print(f'ps {ps} {len(ps)}')
        # Query the region given by the current minimum distance, to see if other points are nearer
        mDist = self.__minDist(ps, p)            
        # Create the region defined by this distance
        ll = Point(p.x - mDist[1], p.y - mDist[1])
        ur = Point(p.x + mDist[1], p.y + mDist[1])            
        ps2, count = self.regionQuery(Bbox([ll, ur]))
        #print(f'ps2 {count} {len(ps)}')                            
        if (count > oldCount): # We only need to recalculate if there are new candidate points
            mDist = self.__minDist(ps2, p)                        
        return mDist
        
    def __minDist(self, points, p):
        dist = []
        for q in points:
            d = p.distEuclidean(q)
            dist.append([q, d]) # Store points and their distances
        mDist = min(dist, key=lambda dist: dist[1]) # Use lamda function to get get minimum distance 
        return mDist

    def bruteNearestPoint(self, p):
        return self.__minDist(self.bigArray, p)

    # Reshape one-dimension list to nCols*nRows list
    def reshapeList(self, list):
        list_new = []
        for i in range(0, self.nRows):
            list_row = []
            for j in range(0, self.nCols):
                list_row.append(list[i*self.nCols+j])
            list_new.append(list_row)
        return list_new

    # Get coordinates of centroids of rasters
    def __getCoords(self):
        coords_x = [self.bBox.ll.x+(self.res/2)]*self.nCols   # List of x coordinates of left bottom raster
        coords_add_x = list(range(0, (self.nCols)*self.res, self.res))
        coords_x = [a+b for a, b in zip(coords_x, coords_add_x)]   # List of x coordinates of one row

        coords_y = [self.bBox.ll.y+(self.res/2)]*self.nRows   # List of y coordinates of left bottom raster
        coords_add_y = list(range(0, (self.nRows)*self.res, self.res))
        coords_y = [a+b for a, b in zip(coords_y, coords_add_y)]   # List of y coordinates of one column

        # Get a nCols*nRows list storing coordinates of rasters centroids
        coords = []
        for i in range(0,self.nRows):
            coords_row = []
            for j in range(0,self.nCols):
                coords_row.append(Point(coords_x[j], coords_y[i]))
            coords.append(coords_row)
        
        # Reshape coords to one-dimension list 
        coords_new = []
        for i in range(0, len(coords)):
            for j in range(0, len(coords[i])):
                coords_new.append(coords[i][j])

        return coords_new
    

                
####### Rasterize Polygon #######
class borderRaster():
    # initialise the index
    def __init__(self, coords, res, border):   # border should be a polygon class
        self.coords = coords   # Coordinates of rasters
        self.res = res   # Resolution
        self.raster_vertex = self.rasterVertex() # Store four vertexes of rasters

        # Get rasters inside border
        self.inside_raster = []  # Store 1 or None, 1 means inside, None means outside
        # Avoid the case when the polygon is small (intersect with rasters but not cover any centroid of raster), check whether the vertexes of rasters are inside the border
        if border.area() < res**2*(len(self.coords)/4):   # If the border area is smaller than quater of rasters areas  
            for i in self.raster_vertex:
                if border.containsPoint(i[0]) == True or border.containsPoint(i[1]) == True or border.containsPoint(i[2]) == True or border.containsPoint(i[3]) == True:
                    self.inside_raster.append(1)
                else:
                    self.inside_raster.append(None)
        # Check whether the centroids of rasters are inside the border
        else:               
            for i in coords:
                if border.containsPoint(i) == True:
                    self.inside_raster.append(1)
                else:
                    self.inside_raster.append(None)
        


    # representation
    def __repr__(self):
        ...
    
    # Get four vertexes of rasters
    def rasterVertex(self):
        raster_vertex = []
        for i in self.coords:
            coords_four = [
                Point(i.x-self.res/2, i.y-self.res/2),
                Point(i.x-self.res/2, i.y+self.res/2),
                Point(i.x+self.res/2, i.y-self.res/2),
                Point(i.x+self.res/2, i.y+self.res/2)
            ]
            raster_vertex.append(coords_four)
        return raster_vertex
    

    # Remove values of rasters outside/inside the border
    def remove_outside(self, list, position):  # position should be 1 or None
        for index, i in enumerate(self.inside_raster):
            if i == position:
                list[index] = None
        return list
    



####### Moran's I #######
class MoranI():
    # initialise
    def __init__(self, data):  # input data should be a nCols*nRows list storing the value in each raster
        self.nCols = len(data[0])
        self.nRows = len(data)
        self.maxIndex = (self.nCols * self.nRows) -1

        # Reshape nCols*nRows list to one-dimension list
        self.data_new = []
        for row in data:
            for num in row:
                self.data_new.append(num)
    
    # representation
    def __repr__(self):
        ...

    # Rook weight matrix
    def __rookWM(self):
        matrix_rook = []
        for i in range(0,self.maxIndex+1):
            matrix_rook_row = []
            for j in range(0,self.maxIndex+1):
                if j == i-1 or j == i+1 or j == i-self.nCols or j == i+self.nCols:
                    matrix_rook_row.append(1)
                else:
                    matrix_rook_row.append(0)
            matrix_rook.append(matrix_rook_row)

        return matrix_rook

    # Queen weight matrix
    def __queenWM(self):
        matrix_queen = []
        for i in range(0,self.maxIndex+1):
            matrix_queen_row = []
            for j in range(0,self.maxIndex+1):
                if j == i-1 or j == i+1 or j == i-self.nCols or j == i+self.nCols or j == i+self.nCols-1 or j == i+self.nCols+1 or j == i-self.nCols-1 or j == i-self.nCols+1:
                    matrix_queen_row.append(1)
                else:
                    matrix_queen_row.append(0)
            matrix_queen.append(matrix_queen_row)

        return matrix_queen

    def globalMoranI(self):
        w = self.__queenWM()
        # Mean of all raster values (except None)
        mean = sum(filter(None, self.data_new))/len(list(filter(None, self.data_new)))
        # Store (x(i) - mean)**2, i is the spatial index of rasters
        data_square = []
        for i in self.data_new:
            if i == None:
                data_square.append(None)
            else:
                data_square.append((i-mean)**2)
        
        upper = 0  # Numerator
        lower = 0  # Denominator

        # Follow the global Moran'I formula
        for i in range(0,len(w)):
            # Sum of weights
            lower += sum(w[i])
            # Get the sum of weighted (x(i) - mean)*(x(j) - mean), i and j are the spatial index of rasters
            for j in range(0,len(w[i])):
                if self.data_new[j] == None or self.data_new[i] == None:
                    upper += 0  # Assume values of rasters that are not counted as 0
                else:
                    upper += w[i][j]*(self.data_new[i]-mean)*(self.data_new[j]-mean)
        # N * sum of weighted (x(i) - mean)*(x(j) - mean), n is the number of rasters
        upper = len(list(filter(None,self.data_new)))*upper
        # weighted sum of (x(i) - mean)**2
        lower = sum(filter(None,data_square))*lower

        moranI_global = upper/lower

        return moranI_global        

    def localMoranI(self):
        w = self.__queenWM()
        # standardize the weight matrix by row
        for i in range(0,len(w)):
            weight_sum = sum(w[i])
            for j in range(0,len(w[i])):
                w[i][j] = w[i][j]/weight_sum
        
        # calculate standardized z scores
        mean = sum(filter(None, self.data_new))/len(list(filter(None, self.data_new)))
        data_square = []
        for i in self.data_new:
            if i == None:
                data_square.append(None)
            else:
                data_square.append((i-mean)**2)
        sd = sqrt(sum(filter(None, data_square))/len(list(filter(None, data_square))))
        z_score = []
        for i in self.data_new:
            if i == None:
                z_score.append(None)
            else:
                z = (i-mean)/sd
                z_score.append(z)

        # Calculate sum of weighted z-score
        wz_sum =  []
        for i in w:
            wz_sum_row = 0   # Store sum of weighted z-score in each row
            for a, b in zip(i, z_score):
                if b == None:
                    wz_sum_row += 0
                else:
                    wz_sum_row += a*b
            wz_sum.append(wz_sum_row)
        
        # Calculate local Moran's I
        moranI_local = []   # Store local Moran's I of each raster
        for a,b in zip(z_score,wz_sum):
            if a == None:
                moranI_local.append(None)
            else:
                moranI_local.append(a*b)

        return moranI_local

    def expectedValue(self):
        # Calculated expected value (-1/(N-1)), N is the number of rasters
        expected = -1/(len(list(filter(None, self.data_new)))-1)
        return expected