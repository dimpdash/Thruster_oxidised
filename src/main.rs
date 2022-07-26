extern crate glutin_window;
extern crate graphics;
extern crate opengl_graphics;
extern crate piston;

use glutin_window::GlutinWindow as Window;
use graphics::Transformed;
use graphics::color::{RED, BLACK, BLUE};
use opengl_graphics::{GlGraphics, OpenGL};
use piston::event_loop::{EventSettings, Events};
use piston::input::{RenderArgs, RenderEvent, UpdateArgs, UpdateEvent};
use piston::window::WindowSettings;
use graphics::context::Context;
use graphics::types::{Color};

use rand::{Rng, ThreadRng};
mod vector;
use self::vector::Vec2;
use std::f64::consts::PI;

pub struct App {
    gl: GlGraphics, // OpenGL drawing backend.
    rotation: f64,  // Rotation for the square.
    fuel: Vec<Particle>,
    bounding_box: Vec<Line>,
    rng: ThreadRng,  
}


#[derive(Debug, Clone)]
pub struct Particle {
    pub color: Color,
    pub pos: Vec2,
    pub mass: f64,
    pub vel: Vec2,
    pub radius: f64,
    pub accel: Vec2,
    pub particle_type: ParticleType,
}

impl Particle {
    pub fn new(color: Color, pos: Vec2, vel: Vec2, mass: f64, particle_type: ParticleType) -> Particle {
        Particle {
            color: color,
            pos: pos,
            vel: vel,
            accel: Vec2::new(0.0, 0.0),
            mass: mass,
            radius: 4.0,
            particle_type: particle_type
        }
    }
    fn render(&self, ctx: Context, gl: &mut GlGraphics) {
        let transform = ctx.transform.trans(0.0,0.0);
        let square = graphics::rectangle::centered_square(self.pos.x, self.pos.y, self.radius);
        graphics::ellipse(self.color, square, transform, gl);
    }
}

pub struct Line {
    pub start: Vec2,
    pub end: Vec2,
    pub color: Color,
    pub width: f64,
    pub line_type: LineType
}

#[derive(Debug, Clone)]
pub enum ParticleType {
    FUEL,
    EXHAUST,
}
pub enum LineType {
    BOUNDING_BOX,
    IGINITER,
}

impl Line {
    pub fn new(start: Vec2, end: Vec2, color: Color, line_type: LineType) -> Line {
        Line {
            start: start,
            end: end,
            color: color, 
            width: 2.0,
            line_type: line_type,
        }
    }

    fn render(&self, ctx: Context, gl: &mut GlGraphics) {
        let transform = ctx.transform.trans(0.0, 0.0);
        graphics::line_from_to(self.color, self.width, [self.start.x, self.start.y], [self.end.x, self.end.y], transform, gl);

    }
}

impl App {
    fn render(&mut self, args: &RenderArgs) {
        use graphics::*;

        self.gl.draw(args.viewport(), |c, gl| {
            // Clear the screen.
            clear(BLACK, gl);

            for f in self.fuel.iter(){
                f.render(c, gl);
            }
            
            for line in self.bounding_box.iter(){
                line.render(c, gl);
            }
            
        });
    }

    fn update(&mut self, args: &UpdateArgs) {
        let delta_t = 0.1;
        //update
        for particle in self.fuel.iter_mut(){
            particle.pos = particle.pos + particle.vel * delta_t;
            particle.vel = particle.vel + particle.accel * delta_t;
        }

       
        //check collision of particles
        for p1_i in 0..self.fuel.len() - 1 {
            let (first, second) = self.fuel.split_at_mut(p1_i+1);
            let p1 = &mut first[p1_i];
            for p2 in second {
                let seperation = p1.pos - p2.pos;
                if seperation.length() < p1.radius + p2.radius {
                    let p1_vel = particle_collision_velocity(p1, p2);
                    p2.vel = particle_collision_velocity(p2, p1);
                    p1.vel = p1_vel;
                }       
            }
        }

        let mut new_particles = vec![];
        //check collision of solid objects
        for particle in self.fuel.iter_mut(){
            for line in self.bounding_box.iter_mut() {
                //from https://monkeyproofsolutions.nl/wordpress/how-to-calculate-the-shortest-distance-between-a-point-and-a-line
                let M = line.end - line.start;

                let t_0 = (particle.pos - line.start).dot(M)/(M.dot(M));
                if 0.0 <= t_0 && t_0 <= 1.0 {
                    let d = (particle.pos - (line.start + t_0 * M)).length();
                    if particle.radius+ particle.vel.length() > d {
                        //Point of collision on line
                        let C2 = line.start + t_0 *M;
                        let sep = particle.pos - C2;
                        particle.vel = particle.vel - 2.0*particle.vel.dot(sep) * (sep)/sep.length().powf(2.0);
                        
                        match (&line.line_type, &particle.particle_type) {
                            (LineType::IGINITER, ParticleType::FUEL) => {
                                    let rand_vel = Vec2 {
                                        x: self.rng.gen_range(0.0, 1.0),
                                        y: self.rng.gen_range(0.0, 1.0)
                                    } * particle.vel.length() * 5.0 ;
                                    
                                    let old_vel = particle.vel;
                                    let old_pos = particle.pos;

                                    let displacement = old_vel.normal()*particle.radius*2.0;

                                    particle.particle_type = ParticleType::EXHAUST;
                                    particle.pos = old_pos + displacement;
                                    particle.radius = particle.radius / 2.0;
                                    particle.color = BLUE;
                                    particle.mass = particle.mass /2.0;
                                    particle.vel = old_vel + rand_vel;

                                    let mut other_particle = particle.clone();
                                    other_particle.vel = old_vel - rand_vel;
                                    other_particle.pos = old_pos - displacement;
                                    new_particles.push(other_particle);
                                }
                            (_,_) => {}
                        }
                    }
                }
            }
        }

        self.fuel.extend( new_particles);

    }
}

fn particle_collision_velocity(p1 : &Particle, p2 : &Particle) -> Vec2 {
    let sep = p1.pos - p2.pos;
    p1.vel - (2.0*p2.mass/(p1.mass + p2.mass))*((p1.vel - p2.vel).dot(sep)as f64)*sep/(sep.length().powf(2.0))
}

fn main() {
    // Change this to OpenGL::V2_1 if not working.
    let opengl = OpenGL::V3_2;

    let window_width = 500.0;
    let window_height = 600.0;

    let initial_window_centre = Vec2 {x : window_width/2.0, y : window_height/2.0};


    // Create a Glutin window.
    let mut window: Window = WindowSettings::new("spinning-square", [window_width, window_height])
        .graphics_api(opengl)
        .exit_on_esc(true)
        .build()
        .unwrap();

    let mut rng = rand::thread_rng();
    
    //create fuel
    let mut fuel = vec![];
    for i in 0..30 {
        for j in 0.. 30 {
            let theta = rng.gen_range(0.0, 2.0 * PI);
            let dir = Vec2::new(theta.cos(), theta.sin());
            let distance: f64 = rng.gen_range(1.0, 150.0);
            let mass = rng.gen_range(100.0, 200.0);
            let momentum = rng.gen_range(20.0,400.0) * distance.sqrt();
            let speed = momentum / mass;
    
            fuel.push(Particle::new(
                RED,
                Vec2 { x: 110.0+(i as f64)*10.0, y: 60.0 + (j as f64)*10.0 },
                dir.normal() * speed,
                10.0,
                ParticleType::FUEL,
            ));
        }

    }

    
    // fuel.push(Particle { 
    //     color: RED, 
    //     pos: Vec2 {x : 480.0, y : 100.0}, 
    //     mass: 100.0, 
    //     vel: Vec2 { x: -1.0, y: 0.0}, 
    //     radius: 40.0,
    //     accel: Vec2 { x: 0.0, y: 0.0 },
    // });

    // fuel.push(Particle { 
    //     color: RED, 
    //     pos: Vec2 {x : 115.0, y : 100.0}, 
    //     mass: 100.0, 
    //     vel: Vec2 { x: -100.0, y: 0.0}, 
    //     radius: 5.0,
    //     accel: Vec2 { x: 0.0, y: 0.0 },
    // });


    let GREY = [0.5,0.5,0.5,1.0];
    
    let mut bounding_box = vec![];
    bounding_box.push(Line::new(
        Vec2 {x:100.0,y:50.0},
        Vec2 {x:100.0,y:400.0},
        GREY,
        LineType::BOUNDING_BOX,
    ));

    bounding_box.push(Line::new(
        Vec2 {x:100.0,y:400.0},
        Vec2 {x:150.0,y:400.0},
        GREY,
        LineType::BOUNDING_BOX,
    ));

    bounding_box.push(Line::new(
        Vec2 {x:350.0,y:400.0},
        Vec2 {x:400.0,y:400.0},
        GREY,
        LineType::BOUNDING_BOX,
    ));

    bounding_box.push(Line::new(
        Vec2 {x:310.0,y:400.0},
        Vec2 {x:200.0,y:400.0},
        GREY,
        LineType::BOUNDING_BOX,
    ));

    bounding_box.push(Line::new(
        Vec2 {x:400.0,y:400.0},
        Vec2 {x:400.0,y:50.0},
        GREY,
        LineType::BOUNDING_BOX,
    ));

    bounding_box.push(Line::new(
        Vec2 {x:100.0,y:50.0},
        Vec2 {x:400.0,y:50.0},
        GREY,
        LineType::BOUNDING_BOX,
    ));

    bounding_box.push(Line::new(
        Vec2 {x:290.0,y:400.0},
        Vec2 {x:260.0,y:450.0},
        RED,
        LineType::IGINITER,        
    ));

    bounding_box.push(Line::new(
        Vec2 {x:230.0,y:400.0},
        Vec2 {x:260.0,y:450.0},
        RED,
        LineType::IGINITER,        
    ));

    bounding_box.push(Line::new(
        Vec2 {x:100.0,y:550.0},
        Vec2 {x:150.0,y:400.0},
        GREY,
        LineType::BOUNDING_BOX,
    ));

    bounding_box.push(Line::new(
        Vec2 {x:400.0,y:550.0},
        Vec2 {x:350.0,y:400.0},
        GREY,
        LineType::BOUNDING_BOX,
    ));

    bounding_box.push(Line::new(
        Vec2 {x:360.0,y:400.0},
        Vec2 {x:330.0,y:420.0},
        GREY,
        LineType::BOUNDING_BOX,
    ));

    bounding_box.push(Line::new(
        Vec2 {x:150.0,y:400.0},
        Vec2 {x:180.0,y:420.0},
        GREY,
        LineType::BOUNDING_BOX,
    ));

    // Create a new game and run it.
    let mut app = App {
        gl: GlGraphics::new(opengl),
        rotation: 0.0,
        fuel: fuel,
        bounding_box,
        rng, 
    };

    let mut events = Events::new(EventSettings::new());
    while let Some(e) = events.next(&mut window) {
        if let Some(args) = e.render_args() {
            app.render(&args);
        }

        if let Some(args) = e.update_args() {
            app.update(&args);
        }
    }
}
