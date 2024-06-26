/*

#[derive(Debug, Default, Clone)]
pub struct Offsets {
    pub offset_area: f64,
    pub perimeter: f64,
    pub contour: Vec<(f64, f64)>,
}

#[derive(Default, Debug, Clone)]
pub struct Polygon {
    pub edges: Vec<Edge>,
    pub vertices: HashMap<usize, Vertex>,
    pub offset_margin: f64,
}

*/

use polygon_offsetting::draw_svg::draw_svg_offset;
use polygon_offsetting::Polygon;
use polygon_offsetting::Offset;

fn example_offsetting() -> Result <(), Box<dyn std::error::Error>> {
	// The initial polygon must be closed
	let points: Vec<(f64, f64)> = vec![
		(0., 0.),
		(100., 0.),
		(40., 20.),
		(100., 20.),
		(100., 60.),
		(95., 78.),
		(55., 40.),
		(0., 60.),
		(0.0, 0.0)
	];

	// The size of our margin offset, if this value is egal to 0 no offsetting will be computed
	let offset_size: f64 = 12.25;
	// Tolerance is the arc segments precision in polygon offsetting (rounded corner)
	let tolerance: f64 = 0.1;

	let mut polygon = Polygon::new(&points, offset_size).map_err(|e| { e })?;
	let offset: Offset = polygon.offsetting(tolerance).map_err(|e| { e })?;

	println!("Initial contour length: {:?}", points.len());
	println!("Offset contour length: {:?}", offset.contour.len());
	println!("Offset area: {:?}", offset.area);
	println!("Offset perimeter: {:?}", offset.perimeter);
    
    draw_svg_offset(
	    &points,
	    &offset,
	    "/examples/svg/",
        "example_offsetting",
	).map_err(|e| { 
	    print!("Error creating offset svg: {:?}", e);
	    e 
	})?;

    Ok(())
}

fn main() {
	match example_offsetting() {
		Ok(_offset_polygon) => { println!("Offset Polygon computed"); },
		Err(e) => { println!("Error Offsetting: {:?}", e); }
	}
}

