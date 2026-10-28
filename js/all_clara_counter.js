import { organic } from './organic.js';

function main() {
    
    // an open path
    let grid = (new organic())
        .points(9)
        .grid(3, 3)
        .selectNth(4)
        .scale(0.6)
        .selectAll()
        .alignLeft()
        .alignTop()
        .translate(55,115)

    // a closed path
    const cx = grid._points[4].x;       // Center X
    const cy = grid._points[4].y;       // Center Y
    const r = 300;        // Radius
    const numPoints = 64; // Smoothness (more points = smoother circle)
    
    // 2. Generate the perimeter points (0 to 2*PI)
    const circlePoints2 = Array.from({ length: numPoints + 1 }, (_, i) => {
        const angle = (i / numPoints) * 2 * Math.PI;
            return [
            cx + r * Math.cos(angle), // X coordinate
            cy + r * Math.sin(angle)  // Y coordinate
        ];
    });

    const circlePoints = Array.from({ length: numPoints + 1 }, (_, i) => {
        const angle = (i / numPoints) * 2 * Math.PI;
        const x = cx + r * Math.cos(angle);
        const y = cy + r * Math.sin(angle);
        return [
            Math.max(cx - r, Math.min(cx + r, x)), // Clamp X to rectangle bounds
            Math.max(cy - r, Math.min(cy + r, y))  // Clamp Y to rectangle bounds
        ];
    });

    
    let centers = grid.copy().selectAll()
    // centers.drawDots('draw');

    let corridor = (new organic())
    	.points(circlePoints);

    let voronois = centers.copy().selectAll()
        .voronoi(corridor);

    var line = d3.line()
            .curve(/*d3.curveCardinalClosed*/ /*d3.curveCatmullRomClosed*/ d3.curveBasisClosed );

    //  Colors are handeled by css classes counter-0, counter-1, counter-2, counter-3, counter-4, counter-5, counter-6, counter-7, counter-8
    //  8  5  2
    //  7  4  1
    //  6  3  0


    let howFullOn = [
    	[],
        [3],
        [6, 3],
        [6, 3, 0],
        [6, 3, 0, 4],
        [6, 3, 0, 4, 7],
        [6, 3, 0, 4, 7, 1],
        [6, 3, 0, 4, 7, 1, 5],
        [6, 3, 0, 4, 7, 1, 5, 8],
        [6, 3, 0, 4, 7, 1, 5, 8, 2],
    ];
  
    function drawCounter(num = 0, voronois) {
        let full = num % howFullOn.length;
        const svg = document.getElementById('draw'); // Assuming you have an SVG element
        svg.innerHTML = ''; // Clear previous paths
        svg.setAttribute("viewBox", "0 0 590 590");
        svg.classList.add('counter-' + full);
        svg.classList.add('counter');
        svg.classList.remove('counter-' + ((full - 1 + howFullOn.length) % howFullOn.length));
        for (let [index, v] of voronois.entries()) {
            let v2 = v.copy().selectAll().increaseResolution(3);
            v2.selectAll().simplify(15).selectAll();
            v2.removeDuplicates();

            const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
            path.classList.add('index-' + index);
            path.setAttribute("d", line(v2._points.map(function(a) { return [a.x, a.y]; })));
            svg.appendChild(path);
        }
    }
    let full = -1;
    setInterval(() => {
        full = full + 1;
        drawCounter(full, voronois);
    }, 1000);
    //drawCounter(5, voronois);
}
main();


// M100,194.83333333333348C100,190.33333333333348,100,182.66666666666683,100.00000000000006,178.83333333333348C100.00000000000011,175.00000000000014,100.00000000000024,175.0000000000001,100.0000000000003,175.00000000000009C100.00000000000038,175.00000000000009,100.00000000000038,175.00000000000009,100.00000000000044,175.00000000000009C100.0000000000005,175.00000000000009,100.00000000000063,175.00000000000003,103.83333333333404,175C107.66666666666742,175,115.3333333333341,175,126.12892532348695,175C136.92451731363983,175,150.84903462727888,175,161.64462661743175,175C172.4402186075846,175,180.10688527425125,175,183.9402186075846,178.83333333333334C187.77355194091794,182.66666666666666,187.77355194091794,190.33333333333334,187.77355194091794,194.83333333333334C187.77355194091794,199.33333333333334,187.77355194091794,200.66666666666666,187.77355194091794,205.16666666666666C187.77355194091794,209.66666666666666,187.77355194091794,217.33333333333334,183.9402186075846,221.16666666666666C180.10688527425125,225,172.4402186075846,225,161.6446266174316,225C150.84903462727863,225,136.9245173136393,225,126.12892532348633,225C115.33333333333333,225,107.66666666666667,225,103.83333333333333,221.16666666666666C100,217.33333333333334,100,209.66666666666666,100,205.1666666666667C100,200.66666666666674,100,199.33333333333346,100,194.83333333333348
// M187.77355194091794,194.83333333333334C187.77355194091794,190.33333333333334,187.77355194091794,182.66666666666666,191.60688527425125,178.83333333333334C195.4402186075846,175,203.10688527425125,175,213.90247853597,175C224.69807179768875,175,238.6225916544596,175,249.41818491617835,175C260.21377817789704,175,267.8804448445637,175,271.7137781778971,178.83333333333334C275.5471115112304,182.66666666666666,275.5471115112304,190.33333333333334,275.5471115112304,194.83333333333334C275.5471115112304,199.3333333333333,275.5471115112304,200.6666666666666,275.5471115112304,205.1666666666666C275.5471115112304,209.66666666666654,275.5471115112304,217.33333333333323,271.7137781778971,221.1666666666666C267.8804448445637,224.9999999999999,260.21377817789704,224.99999999999991,253.64447844791087,224.99999999999991C247.07517871792467,224.99999999999991,241.6032459246189,224.99999999999997,235.03394619463268,224.99999999999997C228.46464646464645,225,220.7979797979798,225,214.8514996987802,225C208.90501959958058,225,204.67872606784806,225,202.56557930198178,225C200.45243253611554,225,200.45243253611554,225,198.33928577024926,225C196.226139004383,225,191.99984547265046,225,189.8866987067842,221.16666666666666C187.77355194091794,217.33333333333334,187.77355194091794,209.66666666666666,187.77355194091794,205.16666666666666C187.77355194091794,200.66666666666666,187.77355194091794,199.33333333333334,187.77355194091794,194.83333333333334