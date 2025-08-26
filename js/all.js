import { organic } from './organic.js';

function main() {

    // an open path
    let centerline = (new organic())
        .points(100)
        .grid(100, 1)
        .selectLast(30)
        .bendUp(200)
        //.selectLast(60)
        //.bendUp(-300);
    // a closed path
    let corridor = centerline.copy()
        .selectAll()
        .corridor(15);

    centerline.draw('draw');
    corridor.draw('draw');

    // distribute 10 points equidistantly along the centerline
    let centers = centerline.copy().selectAll()
        .spline(10 * 2)
        .remove(function (p, i) {
            return i % 2 == 0;
        });
    centers.drawDots('draw');

    let voronois = centers.copy().selectAll()
        .voronoi(corridor);

    var line = d3.line()
            .curve(d3.curveBasisClosed);
    for (let v of voronois) {
        const svg = document.getElementById('draw'); // Assuming you have an SVG element
        const path = document.createElementNS("http://www.w3.org/2000/svg", "path");

        path.setAttribute("d", line(v._points.map(function(a) { return [a.x, a.y]; })));
        path.setAttribute('fill', v._color);
        path.setAttribute('fill-opacity', 0.5);
        svg.appendChild(path);
    }
}
main();
