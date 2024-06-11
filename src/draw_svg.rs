use crate::Offset;

use std::path::Path;
use std::fs::File;
use std::io::Write;

#[derive(Debug, Default)]
pub struct Segm {
    p1: (f64, f64),
    p2: (f64, f64)
}

pub fn draw_svg_offset(
    initial_contour: &Vec<(f64, f64)>,
    offset: &mut Offset,
    dir_path: String,
    filename: String) -> Result <(), Box<dyn std::error::Error>>
{
    let mut segments: Vec<Segm> = Vec::new();
    let mut segms: Vec<Segm> = Vec::new();

    for i in 1..initial_contour.len() {
        let mut sgmt: Segm = Segm::default();

        sgmt.p1.0 = initial_contour[i - 1].0;
        sgmt.p1.1 = initial_contour[i - 1].1;
        if i == initial_contour.len() { 
            sgmt.p2.0 = initial_contour[0].0;
            sgmt.p2.1 = initial_contour[0].1;
            segms.push(sgmt);
            break ;
        }
        sgmt.p2.0 = initial_contour[i].0;
        sgmt.p2.1 = initial_contour[i].1;
        segms.push(sgmt);
    }

    let contour_offset = offset.contour.clone();
    for i in 1..contour_offset.len() {
        let mut sgmt: Segm = Segm::default();

        sgmt.p1.0 = contour_offset[i - 1].0;
        sgmt.p1.1 = contour_offset[i - 1].1;
        if i == contour_offset.len() { 
            sgmt.p2.0 = contour_offset[0].0;
            sgmt.p2.1 = contour_offset[0].1;
            segms.push(sgmt);
            break ;
        }
        sgmt.p2.0 = contour_offset[i].0;
        sgmt.p2.1 = contour_offset[i].1;
        segms.push(sgmt);
    }

    segments.append(&mut segms);

    let mut x_min1 = 0.;
    let mut y_min1 = 0.;
    for i in 0..segments.len() {
        if i == 0 {
            x_min1 = segments[i].p1.0;
            y_min1 = segments[i].p1.1;
        } else {
            if segments[i].p1.0 < x_min1 { x_min1 = segments[i].p1.0; }
            if segments[i].p2.0 < x_min1 { x_min1 = segments[i].p2.0; }
            if segments[i].p1.1 < y_min1 { y_min1 = segments[i].p1.1; }
            if segments[i].p2.1 < y_min1 { y_min1 = segments[i].p2.1; }
        }
    }
    segments.iter_mut().for_each(|p| {
        p.p1.0 -= x_min1;
        p.p2.0 -= x_min1;
        p.p1.1 -= y_min1;
        p.p2.1 -= y_min1;

    });

    let mut max_x: f64 = std::f64::NEG_INFINITY;
    let mut max_y: f64 = std::f64::NEG_INFINITY;
    segments.iter().for_each(|segment| {
        let p = &segment.p1;
        if max_x < p.0 { max_x = p.0 }
        if max_y < p.1 { max_x = p.1 }

        let p = &segment.p2;
        if max_x < p.0 { max_x = p.0 }
        if max_y < p.1 { max_y = p.1 }
    });

    let resize_width =  max_x / 1800.0;
    let resize_height = max_y / 850.0;
    let resize = if resize_width > resize_height { resize_width } else { resize_height };
    segments.iter_mut().for_each(|s| {
        s.p1.0 = s.p1.0 / resize;
        s.p1.1 = s.p1.1 / resize;
        s.p2.0 = s.p2.0 / resize;
        s.p2.1 = s.p2.1 / resize;
    });

    let mut viewbox = vec![
        std::f64::INFINITY,
        std::f64::INFINITY,  
        std::f64::NEG_INFINITY,
        std::f64::NEG_INFINITY
    ];

    segments.iter().for_each(|segment| {
        let p = &segment.p1;
        if viewbox[0] > p.0 { viewbox[0] = p.0 }
        if viewbox[1] > p.1 { viewbox[1] = p.1 }
        if viewbox[2] < p.0 { viewbox[2] = p.0 }
        if viewbox[3] < p.1 { viewbox[3] = p.1 }

        let p = &segment.p2;
        if viewbox[0] > p.0 { viewbox[0] = p.0 }
        if viewbox[1] > p.1 { viewbox[1] = p.1 }
        if viewbox[2] < p.0 { viewbox[2] = p.0 }
        if viewbox[3] < p.1 { viewbox[3] = p.1 }
    });

    let mut txt = format!(
        "<?xml version='1.0' encoding='UTF-8' standalone='no'?>
            <svg xmlns='http://www.w3.org/2000/svg' viewBox='{} {} {} {}' id='export' style='background-color:white'>\n", 
        viewbox[0] - 50., viewbox[1] - 50. , viewbox[2] + 100., viewbox[3] + 100. ).to_string();

    let mut i: usize = 0;
    segments.iter().for_each(|segment| {
        let p1 = segment.p1;
        let p2 = segment.p2;
        
        let mut path = "M ".to_string();
        
        path = format!("{} {} {} ", path, p1.0, (viewbox[3] - viewbox[1]) - (p1.1 - viewbox[1]));
        path = format!("{} {} {} ", path, p2.0, (viewbox[3] - viewbox[1]) - (p2.1 - viewbox[1]));
        
        
        let color = "#000000";
        let width = 1.;

        txt = format!("{}<path style='fill:none;stroke-width:{};stroke:{};' d='{}' id='{}' />\n", txt, width, color, path, 0);
        i += 1;
    });

    txt = format!("{}</svg>", txt);
    let env_path = std::env::current_dir().unwrap().to_str().unwrap().to_string();
    let dir =  env_path.clone() + &dir_path;

    if !(Path::new(&dir_path).exists()) {
        std::fs::create_dir_all(&dir)?;
    }
    let mut file = File::create(format!("{}.svg", dir.clone() + &filename)).map_err(|e| { 
        print!("Error on creating offset svg: {:?}", dir.clone() + &filename);
        e 
    })?;
    file.write_all(txt.as_bytes()).unwrap();
    
    Ok(())
}
