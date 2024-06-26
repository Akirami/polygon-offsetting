# polygon-offsetting

Small crate, no dependencies, to offset a polygon (only margin).

## Features
  - very fast
  - from a simple Vec of 2D coordinates (f64 tuples), compute its offset contour, its perimeter and area
  - the initial polygon must be closed 
  - you must ensure that your polygon does not intersect itself

## Examples

`cargo run --example example_offsetting --release`

```rust
fn example_offsetting() -> Result <(), Box<dyn std::error::Error>> {
	// the initial polygon must be closed
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

	// the size of our margin, if egal to 0 no offsetting will be computed
	let offset_size: f64 = 12.25;
	// tolerance is the arc segments precision in polygon offsetting (rounded corner)
	let tolerance: f64 = 0.1;

	let mut polygon = Polygon::new(&points, offset_size).map_err(|e| { e })?;
	let offset: Offset = polygon.offsetting(tolerance).map_err(|e| { e })?;

	println!("Initial contour length: {:?}", points.len());
	println!("Offset contour length: {:?}", offset.contour.len());
	println!("offset area: {:?}", offset.area);
	println!("offset perimeter: {:?}", offset.perimeter);
    
    draw_svg_offset(
	    &points,
	    &offset,
	    "/examples/svg/",
        "example_offsetting",
	).map_err(|e| { 
	    print!("Error on creating offset svg: {:?}", e);
	    e 
	})?;

    Ok(())
}
```

<img src="./examples/svg/example_offsetting.svg" width = "400">

## Licence
```
MIT License

Copyright (c) 2024 Akirami

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

