[![](https://jitpack.io/v/micycle1/DistanceField.svg)](https://jitpack.io/#micycle1/DistanceField)


# DistanceField

🚧Repo & code under construction🚧

Visualises distance fields for 2D shapes (in Processing).

In this library a distance field defines, for every point (pixel) in a shape, the euclidean distance between it and a single origin point.
Note this is unlike a *signed distance function*, which defines for every point the shortest distance to the shape's boundary.

### Example
![image](https://user-images.githubusercontent.com/9304234/116440088-06323f00-a848-11eb-80a4-872b7b3b7ef7.png)

### Technique
* Triangulate the shape (with refinement) using the shape's *(x, y)* vertices.
* Model the triangulation as a graph; compute the shortest path distance from a given origin vertex to every other vertex.
* Retriangulate the vertices (now *(x, y, z)*, where z is the corresponding shortest path distance).
* Apply an interpolator to the new triangulation, which produces an interpolated distance value for each pixel. 3 different interpolators are provided (offering a speed-quality tradeoff).