# Experiment in visualization

Mouse-human hibrids contained a user interface design inspired by organic space filling blobs. Here an experiment to conver those designs into a package that can be used for other purposes as well. For example in order to display a colormap.

Generate a series of blobs:

![Series of blobs](https://github.com/HaukeBartsch/organic/raw/main/images/bendUp.png)

```javascript
import { organic } from './organic.js';

let centerline = (new organic())
    .points(100)
    .grid(100, 1)
    .selectLast(30)
    .bendUp(200);

let corridor = centerline.copy()
    .selectAll()
    .corridor(15);

// distribute 10 points equidistantly along the centerline
let centers = centerline.copy().selectAll()
    .spline(10 * 2 + 1)
    .remove(function (p, i) {
        return i % 2 == 0;
    });

centerline.draw('draw');
corridor.draw('draw');
centers.drawDots('draw');

let voronois = centers.copy().selectAll()
    .voronoi(corridor);

var line = d3.line()
        .curve(d3.curveBasisClosed);
for (let v of voronois) {
    const svg = document.getElementById('draw');
    const path = document.createElementNS("http://www.w3.org/2000/svg", "path");

    path.setAttribute("d", line(v._points.map(function(a) { return [a.x, a.y]; })));
    path.setAttribute('fill', v._color);
    path.setAttribute('fill-opacity', 0.5);
    svg.appendChild(path);
}
```

![Straight series of blobs](https://github.com/HaukeBartsch/organic/raw/main/images/straight.png)
