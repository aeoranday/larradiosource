# `Geometry` Configuration
The `Geometry` class is an abstract base class with current concrete classes `Cylinder` and `RectangularPrism`.
For the time being, the use of cylinders and rectangular prisms covers many LArTPC designs.

[Geometry UML](docs/uml/geometry.svg "Geometry UML")

Configuring the concrete classes goes as follows:
```toml
[RectangularPrism]
type = "rectangular_prism"
origin = [0, 0, 0]
size = [1, 1, 1]

[Cylinder]
type = "cylinder"
origin = [0, 0, 0]
radius = 1
height = 1
```

The `Geometry` class is not used on it's own.
More frequently, a user would configure a `Detector` or `Source` geometry using:
```toml
[source_name.geometry]
type = "cylinder"
origin = [0, 0, 0]
radius = 1
height = 1

[detector_name.geometry]
type = "rectangular_prism"
origin = [0, 0, 0]
size = [1, 1, 1]
```

For both volumes, the extensions (`radius` and `height` or `size`) start from the `origin`.
In the above examples, the rectangular prism would be a unit cube with one corner at `(0, 0, 0)` and another corner at `(1, 1, 1)` and the cylinders would be a unit circle centered at `(0, 0, 0)` and extending by 1 unit.
