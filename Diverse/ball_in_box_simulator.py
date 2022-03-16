# SOURCES :
#

import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import turtle

def make_rectangle(x0,y0,W,H,color="green"):
    walls = turtle.Turtle()
    turtle.speed(0)
    walls.color(color)
    walls.penup()
    walls.setx(x0-W/2)
    walls.sety(y0+H/2)
    walls.pendown()
    walls.goto(x0+W/2,x0+H/2)
    walls.goto(x0+W/2,x0-H/2)
    walls.goto(x0-W/2,x0-H/2)
    walls.goto(x0-W/2,x0+H/2)
    walls.begin_poly()
    walls.hideturtle()
    return walls;

def make_ball(d,speed=0,x0=0,y0=0,color="white"):

    ball = turtle.Turtle()
    ball.setx(x0)
    ball.sety(y0)
    ball.shape("circle")
    ball.color(color)
    ball.shapesize(d)  # set radius
    ball.penup()
    ball.speed(speed)

    return ball;

def gravitation_effect(v,dt,g=9.81):
    v[1] -= dt*g;

#make screen
wn = turtle.Screen()
wn.bgcolor("black")
wn.title("Ball in Box")
wn.tracer(0) #dont update by itself

#make ball
ball1 = make_ball(d=0.5,speed=1,x0=0,y0=0,color="green")
#box
W = 600
H = 500
x0 = 0
y0 = 0
wall = make_rectangle(0,0,W,H,"red")

#physical frame
v = [1,1] #m/s
sf = 1
dt = 0.001
g = 9.81
friction = 0.1 # slows speed
t = 0

#turtle.mainloop()
ball1.pendown() #draw trajectory
while True:
    wn.update() #update canvas
    ball1.sety(ball1.ycor()+v[1]*sf)
    ball1.setx(ball1.xcor()+v[0]*sf)

    if(ball1.xcor() <= -(x0+W/2) or ball1.xcor() >= (x0+W/2)): #leftbounce
        v[0] *= -1
        v = np.multiply(v,1-friction)

    if(ball1.ycor() <= -(y0+H/2) or ball1.ycor() >= (y0+H/2)):
        v[1] *= -1
        v = np.multiply(v,1-friction)

    gravitation_effect(v,dt)
