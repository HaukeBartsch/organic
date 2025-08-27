import SplayTree from './splay.esm.js';
import { default as polygonClipping} from './polygon-clipping.esm.js';

// Code to solve a system of linear equations using Gauss-Jordan elimination.
// The main entry point is the solve() function.
//
// Eli Bendersky [https://eli.thegreenplace.net]
// This code is in the public domain.
'use strict';

// This code uses an array-of-arrays representation of 2D matrices, e.g.:
// 
// let mat = [
//     [-1, 4, -2, -15],
//     [-4, 6, 1, -5],
//     [-6, -6, -2, -10],
// ];

// solve solves the system of linear equations Ax = b, where A is a matrix
// and b is an array representing a column vector. The solution x is returned
// as an array. solve throws an exception if the system doesn't have a unique
// solution.
// A is modified in place - it should be cloned outside this function if you
// want to preserve the original.
function solve(A, b) {
    // Step 1: create the augmented matrix [A|b], while making sure all
    // dimensions match. The resulting matrix has R rows and R+1 columns.
    let R = A.length;
    if (R != b.length) {
        throw new Error("A and b must have the same number of rows");
    }
    for (let i = 0; i < R; i++) {
        if (A[i].length != R) {
            throw new Error("A must be square");
        }
        A[i].push(b[i]);
    }

    // Step 2: perform Gaussian elimination on the augmented matrix. This
    // modifies A to be in row echelon form.
    gaussEliminate(A);

    // Step 3: back-substitution. This modifies A to be in reduced row
    // echelon form (Gauss-Jordan elimination).
    for (let i = R - 1; i >= 0; i--) {
        // For each row, take its pivot and divide the last column by it,
        // then eliminate the pivot from all rows above.
        let pivot = A[i][i];
        if (pivot == 0) {
            throw new Error("System has no unique solution");
        }
        for (let j = i - 1; j >= 0; j--) {
            let f = A[j][i] / pivot;
            A[j][i] = 0;
            A[j][R] -= A[i][R] * f;
        }
        A[i][i] = 1;
        A[i][R] /= pivot;
    }

    // Step 4: extract the solution vector from the last column of A.
    let x = [];
    for (let i = 0; i < R; i++) {
        x.push(A[i][R]);
    }
    return x;
}

// findPivotRow finds the "pivot" row in arr, for column col and beginning
// with startRow. The pivot row is the row with the largest (in absolute value)
// element in column col among rows [startRow:arr.length). The index of the
// pivot row is returned.
function findPivotRow(arr, startRow, col) {
    let maxidx = startRow;
    for (let i = startRow + 1; i < arr.length; i++) {
        if (Math.abs(arr[i][col]) > Math.abs(arr[maxidx][col])) {
            maxidx = i;
        }
    }
    return maxidx;
}

// swapRows swaps rows i and j in arr, in place.
function swapRows(arr, i, j) {
    if (i != j) {
        let tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
}

// gaussEliminate performs Gaussian elimination on arr, in place. After running,
// arr will be in row echelon form. It operates on arbitraryly sized matrices.
// This code follows the pseudocode from Wikipedia, with partial pivoting
// (https://en.wikipedia.org/wiki/Gaussian_elimination). It selects the largest
// possible absolute value for each column to improve numerical stability.
function gaussEliminate(arr) {
    let nrows = arr.length;
    let ncols = arr[0].length;

    let h = 0, k = 0;

    while (h < nrows && k < ncols) {
        // Find the pivot row for column k.
        let pivotRow = findPivotRow(arr, h, k);
        if (arr[pivotRow][k] == 0) {
            // No pivot in this column; move on to the next one.
            k++;
        } else {
            // Swap current row with the pivot row, so we can use the pivot's
            // leading element to eliminate below.
            swapRows(arr, h, pivotRow);

            for (let i = h + 1; i < nrows; i++) {
                let f = arr[i][k] / arr[h][k];
                arr[i][k] = 0;
                for (let j = k + 1; j < ncols; j++) {
                    arr[i][j] -= arr[h][j] * f;
                }
            }
            h++;
            k++;
        }
    }
}

function buildSplineEquations(xs, ys) {
    // Npolys is the number of (cubic) polynomials we interpolate between
    // the given points. Ncoeffs is the number of coefficients they all have
    // together (4 per poly: ax^3 + bx^2 + cx + d).
    // Npoints is the number of points.
    const Npoints = xs.length;
    const Npolys = Npoints - 1;
    const Ncoeffs = 4 * Npolys;

    // The matrix A is the coefficient matrix for the system of linear
    // equations we need to solve to find the coefficients of the polynomials.
    // It has Ncoeffs rows and columns.
    // for each poly i, A[4*i][0..3] are the 4 coefficients of this poly.
    let A = [];
    for (let i = 0; i < Ncoeffs; i++) {
        A.push(Array(Ncoeffs).fill(0));
    }

    // The vector b is the right-hand side of the system of linear equations.
    // It has Ncoeffs values.
    let b = Array(Ncoeffs).fill(0);

    // Now we start filling in the matrix A and vector b.
    // First, we fill in the constraints that the polynomials must pass
    // through the given points. This populates the first 2*Npolys rows.
    let nrow = 0;
    for (let i = 0; i < Npolys; i++, nrow += 2) {
        // Poly i passes through points i and i+1.
        A[nrow][4 * i] = xs[i] ** 3;
        A[nrow][4 * i + 1] = xs[i] ** 2;
        A[nrow][4 * i + 2] = xs[i];
        A[nrow][4 * i + 3] = 1;
        b[nrow] = ys[i];

        A[nrow + 1][4 * i] = xs[i + 1] ** 3;
        A[nrow + 1][4 * i + 1] = xs[i + 1] ** 2;
        A[nrow + 1][4 * i + 2] = xs[i + 1];
        A[nrow + 1][4 * i + 3] = 1;
        b[nrow + 1] = ys[i + 1];
    }

    // Constraints for the first derivatives. This works on non-boundary points,
    // so it gives us (Npolys - 1) equations.
    for (let i = 0; i < Npolys - 1; i++, nrow++) {
        // Poly i and poly i+1 must have the same first derivative at
        // point i+1.
        A[nrow][4 * i] = 3 * xs[i + 1] ** 2;
        A[nrow][4 * i + 1] = 2 * xs[i + 1];
        A[nrow][4 * i + 2] = 1;
        A[nrow][4 * (i + 1)] = -3 * xs[i + 1] ** 2;
        A[nrow][4 * (i + 1) + 1] = -2 * xs[i + 1];
        A[nrow][4 * (i + 1) + 2] = -1;
    }

    // Constraints for the second derivatives. This also gives us (Npolys - 1)
    // equations.
    for (let i = 0; i < Npolys - 1; i++, nrow++) {
        // Poly i and poly i+1 must have the same second derivative at
        // point i+1.
        A[nrow][4 * i] = 6 * xs[i + 1];
        A[nrow][4 * i + 1] = 2;
        A[nrow][4 * (i + 1)] = -6 * xs[i + 1];
        A[nrow][4 * (i + 1) + 1] = -2;
    }

    // The final two equations come from the "natural" boundary conditions;
    // the first and last polys must have zero second derivative at the
    // endpoints.
    A[nrow][0] = 6 * xs[0];
    A[nrow][1] = 2;
    A[nrow + 1][4 * (Npolys - 1)] = 6 * xs[Npolys];
    A[nrow + 1][4 * (Npolys - 1) + 1] = 2;

    return [A, b];
}

// linspace returns an array of numPoints values distributed linearly in
// the (inclusive) range [start,end], just like Numpy's linspace.
function linspace(start, end, numPoints) {
    if (numPoints === undefined || numPoints < 2) {
        return [start, end];
    }

    const step = (end - start) / (numPoints - 1);
    return new Array(numPoints).fill(null).map((_, i) => start + i * step);
}

// calcSinc calculates the sinc(x) function and returns a y value.
// https://en.wikipedia.org/wiki/Sinc_function
function calcSinc(x) {
    if (x == 0) {
        return 1
    } else {
        return Math.sin(Math.PI * x) / (Math.PI * x);
    }
}


function doInterpolate(xs, ys, N) {
    // Perform interpolation on xs, ys to get the coefficients of the splines.
    let [A, b] = buildSplineEquations(xs, ys);
    let coeffs = solve(A, b);
    console.log(coeffs);

    // Create N points linearly spaced between the min and max of xs, and
    // calculate the corresponding py for each px using the appropriate curve.
    let pxs = linspace(Math.min(...xs), Math.max(...xs), N);

    let pys = Array(N).fill(0);
    for (let i = 0; i < N; i++) {
        let px = pxs[i];
        // Find the number of the curve for px, based on which points from
        // xs it's between. Can be done more efficiently with binary
        // search, but this is good enough for a demo.
        let curveIndex = -1;
        for (let j = 0; j < xs.length - 1; j++) {
            // is px between xs[j] and xs[j+1]? If yes, we found the curve!
            if (px >= xs[j] && px <= xs[j + 1]) {
                curveIndex = j;
                break;
            }
        }
        if (curveIndex < 0) {
            alert(`curve index not found for xs[${i}]=${xs[i]}`);
        }

        // With the curve index in hand, we can calculate py based on the
        // relevant curve coefficients from coeffs.
        let [a, b, c, d] = coeffs.slice(curveIndex * 4, curveIndex * 4 + 4);
        pys[i] = a * px ** 3 + b * px ** 2 + c * px + d;
    }

    return [pxs, pys];
}

function createSvgPath(points, closePath = false) {
    if (points.length === 0) {
        return "";
    }

    let pathData = `M ${points[0].x} ${points[0].y}`; // Move to the first point

    for (let i = 1; i < points.length; i++) {
        pathData += ` L ${points[i].x} ${points[i].y}`; // Draw lines to subsequent points
    }

    if (closePath) {
        pathData += " Z"; // Close the path
    }

    return pathData;
}


class organic {
    constructor() {
        this._points = [];
        this._selected = [];
        this._closed = false;
        this.setColor(-1);
    }

    setColor(idx = -1) {
        if (idx < 0)
            idx = Math.floor(Math.random() * 100);
        let ar = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6"];
        this._color = ar[idx % ar.length];
        return this;
    }

    points(numPoints) {
        for (let i = 0; i < numPoints; i++) {
            this._points.push({ x: 0, y: 0 });
        }
        return this;
    }

    remove(index) {
        if (typeof index === 'function') {
            this._points = this._points.filter((p, i) => !index(p, i));
            return this;
        }

        if (index < 0)
            index = this._points.length + index;
        if (index >= 0 && index < this._points.length) {
            this._points.splice(index, 1);
        }
        return this;
    }

    copy() {
        let newObj = new organic();
        newObj._points = JSON.parse(JSON.stringify(this._points));
        newObj._selected = JSON.parse(JSON.stringify(this._selected));
        newObj._closed = this._closed;
        newObj.setColor(); // new random color   
        return newObj;
    }

    draw(elementId) {
        const svg = document.getElementById(elementId); // Assuming you have an SVG element
        const path = document.createElementNS("http://www.w3.org/2000/svg", "path");

        const dAttribute = createSvgPath(this._points, this._closed); // Create a closed path
        path.setAttribute("d", dAttribute);
        path.setAttribute("stroke", this._color);
        path.setAttribute("fill", "none"); // Or a color if you want to fill it
        svg.appendChild(path);
        return path;
    }

    drawDots(elementId, radius = 5) {
        const svg = document.getElementById(elementId); // Assuming you have an SVG element
        const group = document.createElementNS("http://www.w3.org/2000/svg", "g");

        for (let p of this._points) {
            if (this._selected.length > 0 && !this._selected.includes(this._points.indexOf(p)))
                continue;
            var circle = document.createElementNS("http://www.w3.org/2000/svg", 'circle');
            circle.setAttributeNS(null, 'cx', p.x);
            circle.setAttributeNS(null, 'cy', p.y);
            circle.setAttributeNS(null, 'r', radius);
            circle.setAttributeNS(null, 'fill', this._color);
            group.appendChild(circle);
        }
        svg.appendChild(group);
        return group;
    }

    drawIds(elementId, fontsize = 5) {
        const svg = document.getElementById(elementId); // Assuming you have an SVG element
        const group = document.createElementNS("http://www.w3.org/2000/svg", "g");

        for (let p of this._points) {
            if (this._selected.length > 0 && !this._selected.includes(this._points.indexOf(p)))
                continue;
            var text = document.createElementNS("http://www.w3.org/2000/svg", 'text');
            text.setAttributeNS(null, 'x', p.x);
            text.setAttributeNS(null, 'y', p.y);
            text.setAttributeNS(null, 'dx', fontsize/2);
            text.setAttributeNS(null, 'dy', -fontsize/2);
            text.setAttributeNS(null, 'font-size', fontsize);
            // circle.setAttributeNS(null, 'r', radius);
            text.setAttributeNS(null, 'fill', this._color);
            var textNode = document.createTextNode(this._points.indexOf(p));
            text.appendChild(textNode);
            group.appendChild(text);
        }
        svg.appendChild(group);
        return group;
    }


    grid(cols, rows) {
        // assume the screen is 1000x1000 
        let screenWidth = 1000;
        let screenHeight = 1000;

        let arrangedPoints = [];
        let spacingX = 800 / (cols - 1);
        let spacingY = 600 / (rows - 1);
        if (rows < 2)
            spacingY = 0;
        if (cols < 2)
            spacingX = 0;
        let shiftX = (screenWidth / 2) - (800 / 2);
        let shiftY = (screenHeight / 2) - (600 / 2);

        for (let i = 0; i < cols; i++) {
            for (let j = 0; j < rows; j++) {
                arrangedPoints.push({ x: i * spacingX + shiftX, y: j * spacingY + shiftY });
            }
        }
        this._points = arrangedPoints;
        return this;
    }

    selectLast(n) {
        this._selected = [];
        for (let i = Math.max(0, this._points.length - n); i < this._points.length; i++) {
            this._selected.push(i);
        }
        return this;
    }

    selectFirst(n) {
        this._selected = [];
        for (let i = 0; i < Math.min(n, this._points.length); i++) {
            this._selected.push(i);
        }
        return this;
    }

    bendUp(amount = 1) {
        let counter = 0;
        let scale = -amount;
        for (let i of this._selected) {
            let temp = this._points[i];
            this._points[i] = { x: temp.x, y: temp.y + (Math.pow((counter / this._selected.length), 2.0)) * scale };
            counter++;
        }
        return this;
    }

    selectAll() {
        this._selected = [];
        for (let i = 0; i < this._points.length; i++) {
            this._selected.push(i);
        }
        return this;
    }

    corridor(radius = 10) {
        // use the line normal to compute offset points, connect those to a connected line
        let newPointsUp = [];
        let newPointsDown = [];
        if (this._selected.length == 0) {
            this.selectAll();
        }
        for (let i = 0; i < this._selected.length; i++) {
            let p0 = this._points[i];
            let p1 = null;
            if (i == this._selected.length - 1) // continue the line
                p1 = {
                    x: this._points[i].x + (this._points[i].x - this._points[i - 1].x),
                    y: this._points[i].y + (this._points[i].y - this._points[i - 1].y)
                };
            else
                p1 = this._points[i + 1];
            let normal = { x: p1.y - p0.y, y: -(p1.x - p0.x) };
            let length = Math.sqrt(normal.x * normal.x + normal.y * normal.y);
            normal.x /= length;
            normal.y /= length;

            newPointsUp.push({ x: p0.x - (normal.x * radius), y: p0.y - (normal.y * radius) });
            newPointsDown.push({ x: p0.x + (normal.x * radius), y: p0.y + (normal.y * radius) });
        }
        // we should add one more point... 
        this._points = newPointsUp.concat(newPointsDown.reverse());
        this._closed = true;
        return this;
    }

    spline(resolution) {
        const lineGenerator = d3.line()
        .x(d => d.x)
        .y(d => d.y)
        .curve(d3.curveBasis);
        let svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
        let path = document.createElementNS("http://www.w3.org/2000/svg", "path");
        path.setAttribute("d", lineGenerator(this._points));
        let curveLength = path.getTotalLength();
        let newPoints = [];
        for (let i = 0; i < resolution; i++) {
            let p = path.getPointAtLength( (i / resolution) * curveLength);
            newPoints.push({ x: p.x, y: p.y });
        }
        let p = path.getPointAtLength(1.0 * curveLength);
        newPoints.push({ x: p.x, y: p.y });
        this._points = newPoints;
        return this;
    }

    removeDuplicates(tolerance = 1e-4) {
        let newPoints = [];
        let lastPoint = null;
        for (let p of this._points) {
            if (lastPoint == null) {
                newPoints.push(p);
                lastPoint = p;
                continue;
            }
            let dx = p.x - lastPoint.x;
            let dy = p.y - lastPoint.y;
            let dist = Math.sqrt(dx * dx + dy * dy);
            if (dist > tolerance) {
                newPoints.push(p);
                lastPoint = p;
            }
        }
        this._points = newPoints;
        if (this._closed) {
            // the first and last point also should be different, remove last point until they are
            while (this._points.length > 1) {
                let p0 = this._points[0];
                let p1 = this._points[this._points.length - 1];
                let dx = p0.x - p1.x;
                let dy = p0.y - p1.y;
                let dist = Math.sqrt(dx * dx + dy * dy);
                if (dist > tolerance)
                    break;
                this._points.pop();
            }
        }

        this.selectAll();
        return this;
    }

    simplify(tolerance = 1.0) {
        let numPointsBefore = 0;
        do {
            numPointsBefore = this._points.length;
            this.selectAll()._simplify(tolerance).selectAll();
        } while (this._points.length < numPointsBefore);
        return this;
    }

    // we call this several times until no more change
    _simplify(tolerance = 1.0) {
        let pointsToSimplify = [];
        let smallestDistance = Infinity;
        for (let i of this._selected) {
            pointsToSimplify.push([this._points[i].x, this._points[i].y]);
            if (pointsToSimplify.length > 1) {
                let dx = pointsToSimplify[pointsToSimplify.length - 1][0] - pointsToSimplify[pointsToSimplify.length - 2][0];
                let dy = pointsToSimplify[pointsToSimplify.length - 1][1] - pointsToSimplify[pointsToSimplify.length - 2][1];
                let dist = Math.sqrt(dx * dx + dy * dy);
                if (dist < smallestDistance)
                    smallestDistance = dist;
            }
        }
        if (tolerance < smallestDistance)
            return this;

        let newPoints = [];
        // todo: remove points that are close to both successor and predecessor
        newPoints.push(pointsToSimplify[0]);
        for (let i = 1; i < pointsToSimplify.length-1; i++) {
            let dx = pointsToSimplify[i][0] - pointsToSimplify[i+1][0];
            let dy = pointsToSimplify[i][1] - pointsToSimplify[i+1][1];
            let dist = Math.sqrt(dx * dx + dy * dy);
            let dx2 = pointsToSimplify[i][0] - pointsToSimplify[i-1][0];
            let dy2 = pointsToSimplify[i][1] - pointsToSimplify[i-1][1];
            let dist2 = Math.sqrt(dx2 * dx2 + dy2 * dy2);
            if ((dist + dist2) < tolerance) { // skip the current point
                newPoints.push(pointsToSimplify[i+1]);
                i++;
                continue; // skip this point
            }
            newPoints.push(pointsToSimplify[i]);
        }
        newPoints.push(pointsToSimplify[pointsToSimplify.length - 1]); // always keep the last point
        this._points = newPoints.map(p => { return { x: p[0], y: p[1] }; });
        this.selectAll();
        return this;
    }

    spline2(resolution) { // this only works for graphs, not if the splines curves in on itself
        let origXs = [];
        let origYs = [];
        for (let i of this._selected) {
            let p = this._points[i];
            origXs.push(p.x);
            origYs.push(p.y);
        }
        let [pxs, pys] = doInterpolate(origXs, origYs, resolution);
        this._points = [];
        for (let i = 0; i < pxs.length; i++) {
            this._points.push({ x: pxs[i], y: pys[i] });
        }
        // Placeholder for spline logic
        return this;
    }

    rotateIndex(steps = 1) {
        if (steps == 0)
            return this;
        function rotateArray(arr, k, direction = 'left') {
            const n = arr.length;
            k = k % n;

            if (k === 0) return [...arr]; // Return a shallow copy if no rotation

            if (direction === 'left') {
                return arr.slice(k).concat(arr.slice(0, k));
            } else if (direction === 'right') {
                return arr.slice(n - k).concat(arr.slice(0, n - k));
            } else {
                throw new Error("Direction must be 'left' or 'right'.");
            }
        }
        if (steps > 0)
            this._points = rotateArray(this._points, steps, 'left');
        else
            this._points = rotateArray(this._points, -steps, 'right');
        return this;
    }

    // add a new point between each existing point
    increaseResolution(factor = 1) {
        let newPoints = [];
        for (let i = 0; i < this._points.length - 1; i++) {
            let p0 = this._points[i];
            let p1 = this._points[i + 1];
            newPoints.push(p0);
            if (i in this._selected) {
                for (let j = 1; j <= factor; j++) {
                    let t = j / (factor + 1);
                    newPoints.push({ x: p0.x * (1 - t) + p1.x * t, y: p0.y * (1 - t) + p1.y * t });
                }
            }
        }
        newPoints.push(this._points[this._points.length - 1]);
        this._points = newPoints;
        return this;
    }

    // border is the polygon used to clip
    voronoi(border) {
        //var polygons = d3.geom.voronoi(this._points);

        var isCornerRadiusAbsolute = true;
        var cornerRadius = 33;
        function resampleSegments(points) {
            if (points.length == 0)
                return points;

            let i = -1;
            let n = points.length;
            let p0 = points[n - 1];
            let x0 = p0[0];
            let y0 = p0[1];
            let p1, x1, y1;
            let points2 = [];

            while (++i < n) {
                p1 = points[i];
                x1 = p1[0];
                y1 = p1[1];

                let finalRadius = 0;

                if (isCornerRadiusAbsolute) {
                    let distance = Math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2);
                    let distFromPoint = cornerRadius / distance;
                    finalRadius = distFromPoint >= 0.5 ? 0.5 : distFromPoint;
                }
                else {
                    finalRadius = cornerRadius;
                }

                points2.push(
                    [x0 + (x1 - x0) * finalRadius,
                    y0 + (y1 - y0) * finalRadius],
                    [x0 + (x1 - x0) * (1 - finalRadius),
                    y0 + (y1 - y0) * (1 - finalRadius)],
                    p1
                );

                p0 = p1;
                x0 = x1;
                y0 = y1;
            }
            return points2;
        }

        // wrong orientation, the reverse fixes the problem
        function circle(cx, cy, r, n) {
            var points = [];
            d3.range(1e-6, 2 * Math.PI, 2 * Math.PI / n).map(function (θ, i) {
                var point = [cx + Math.cos(θ) * r, cy + Math.sin(θ) * r];
                points.push(point);
                return point;
            });
            points = points.reverse();
            return points;
        }
        //var line = d3.line()
        //    .curve(d3.curveBasisClosed)

        var bounds = d3.geom.polygon(border._points.map(p => [p.x, p.y]));

        var weightedVoronoi = d3.weightedVoronoi()
            .x(function (d) { return d[0]; })                     // set the x coordinate accessor
            .y(function (d) { return d[1]; })                     // set the y coordinate accessor
            .weight(function (d) {
                return 1;
            })  // set the weight accessor
            .clip(bounds);  // set the clipping polygon

        //var ps = d3.selectAll('.point').data(this._points.map(p => [p.x, p.y]));    
        var vor = weightedVoronoi(this._points.map(p => [p.x, p.y]));                        // compute the weighted Voronoi tessellation
        var newPolygons = [];
        for (let cell of vor) {
            if (cell.length > 2) {
                // this clipping only uses the convex approximation of the border, do better using polygon-clipping
                var erg = polygonClipping.intersection([border._points.map(p => [p.x, p.y])], [cell]);

                // let resampled = resampleSegments(erg[0][0].filter((p, i) => i < erg[0][0].length - 1));
                let resampled = erg[0][0]; // resampleSegments(erg[0][0]);
                let newPoly = new organic();
                for (let p of resampled) {
                    newPoly._points.push({ x: p[0], y: p[1] });
                }
                newPoly._closed = true;
                //newPoly.setColor(newPolygons.length);
                newPolygons.push(newPoly);
            }
        }

        // sort the polygons based on the original points in this._points (closest distance to center of mass is fine)
        let center_of_mass = newPolygons.map(polygon => {
            let comx = 0;
            let comy = 0;
            for (let p of polygon._points) {
                comx += p.x;
                comy += p.y;
            }
            comx /= polygon._points.length;
            comy /= polygon._points.length;
            return { x: comx, y: comy };
        });

        let newPolygons2 = [];
        for (let i = 0; i < this._points.length; i++) {
            let closestIndex = -1;
            let closestDistance = Infinity;
            for (let j = 0; j < center_of_mass.length; j++) {
                let dx = this._points[i].x - center_of_mass[j].x;
                let dy = this._points[i].y - center_of_mass[j].y;
                let dist = Math.sqrt(dx * dx + dy * dy);
                if (dist < closestDistance) {
                    closestDistance = dist;
                    closestIndex = j;
                }
            }
            if (closestIndex >= 0) {
                newPolygons2.push(newPolygons[closestIndex]);
                newPolygons2[newPolygons2.length-1].setColor(newPolygons2.length);
            }
        }

        // returns an array of organic objects
        return newPolygons2;
    }

}
export { organic };