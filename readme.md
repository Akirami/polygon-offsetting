# polygon-offsetting

Small crate, no dependencies, to offset a polygon (only margin).

## Features
  - from a simple Vec of 2D coordinates (f64 tuples), return its offset contour, perimeter and area.
  - the initial polygon must be closed 
  - you must ensure your polygon is not self intersecting

## Examples

`cargo run --example example_offsetting --release`

```rust
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
	let offset_size: f64 = 2.25;
	// Tolerance is the arc segments precision in polygon offsetting (rounded corner)
	let tolerance: f64 = 0.1;

	let mut polygon = Polygon::new(&points, offset_size).map_err(|e| { e })?;
	let offset: Offsets = polygon.offsetting(tolerance).map_err(|e| { e })?;

	println!("initial contour length: {:?}", points.len());
	println!("Offset contour length: {:?}", offset.contour.len());
	println!("offset area: {:?}", offset.offset_area);
	println!("offset perimeter: {:?}", offset.perimeter);
    
    draw_svg_offset(
	    &points,
	    &mut offset.clone(),
	    "/examples/svg/".to_string(),
        "example_offsetting".to_string(),
	).map_err(|e| { 
	    print!("Error on creating offset svg: {:?}", e);
	    e 
	})?;

    Ok(())
}
```

![Example svg output](./examples/svg/example_offsetting.svg)

