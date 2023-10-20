from vpython import *
import random
from scipy.spatial import distance
from sklearn.metrics import mean_squared_error
import math
import matplotlib.pyplot as plt
import numpy as np
import time

AVOGADRO_CONSTANT = 6.02 * math.pow(10,23)
BOLTZMANN_CONSTANT = 1.38 * math.pow(10,-23)

PARTICLE_ATOMIC_MASS = 0.222#kg/mol -- #WE ARE ASSUMING THE PARTICLES ARE RADON
PARTICLE_DIAMETRE = 5 * math.pow(10,-10) #m
PARTICLE_MASS = PARTICLE_ATOMIC_MASS / AVOGADRO_CONSTANT #kg

NUMBER_OF_PARTICLES = int()
VOLUME = float()#m^3
TEMPERATURE = float()  #degrees kelvin

PARTICLE_SPEED = float() #m * s^-1
SIDE_LENGTH = float()#m
TOTAL_WALL_AREA = float()#m^2
WALL_WIDTH = 2 #arbitrary

FRAME_RATE = 60#fs^-1
MAKE_TRAILS = True

particles = list()

left_wall = box()
right_wall = box()
top_wall = box()
bottom_wall = box()
back_wall = box()

global_particle_velocity_change = float() #total count of all changes in velocity, as a result of colliding with the walls



def reflect(particle, axis):
    global global_particle_velocity_change
    if axis == 'x':
        global_particle_velocity_change += (abs(2 * particle.vector.x))
        particle.vector.x = -particle.vector.x
    elif axis == 'y':
        global_particle_velocity_change += (abs(2 * particle.vector.y))
        particle.vector.y = -particle.vector.y
    elif axis == 'z':
        global_particle_velocity_change += (abs(2 * particle.vector.z))
        particle.vector.z = -particle.vector.z

def spawn_particles():
    particles.clear()
    for i in range(NUMBER_OF_PARTICLES):
        particle_instance = Particle(random.uniform(-20,20),random.uniform(-20,20),random.uniform(-20,20),PARTICLE_SPEED)
        particles.append(particle_instance)

def spawn_walls():
    global left_wall,right_wall,top_wall,top_wall,bottom_wall,back_wall
    wall_dict = {
        left_wall : [vector(-SIDE_LENGTH/2-(WALL_WIDTH/2),0,0), vector(WALL_WIDTH,SIDE_LENGTH,SIDE_LENGTH)],
        right_wall : [vector(SIDE_LENGTH/2+(WALL_WIDTH/2),0,0), vector(WALL_WIDTH,SIDE_LENGTH,SIDE_LENGTH)],
        top_wall : [vector(0,SIDE_LENGTH/2+(WALL_WIDTH/2),0),vector(SIDE_LENGTH,WALL_WIDTH,SIDE_LENGTH) ],
        bottom_wall : [vector(0,-SIDE_LENGTH/2-(WALL_WIDTH/2),0),vector(SIDE_LENGTH,WALL_WIDTH,SIDE_LENGTH)],
        back_wall : [vector(0,0,-SIDE_LENGTH/2-(WALL_WIDTH/2)),vector(SIDE_LENGTH,SIDE_LENGTH,WALL_WIDTH)]
    }
    for wall, position_info in wall_dict.items():
        wall.pos = position_info[0]
        wall.size = position_info[1]

class Particle(sphere):
    def __init__(self, x_component,y_component,z_component,_speed):
        super().__init__(pos=vector(random.uniform((-SIDE_LENGTH/2)+PARTICLE_DIAMETRE,(SIDE_LENGTH/2)-PARTICLE_DIAMETRE),random.uniform((-SIDE_LENGTH/2)+PARTICLE_DIAMETRE,(SIDE_LENGTH/2)-PARTICLE_DIAMETRE),random.uniform((-SIDE_LENGTH/2)+PARTICLE_DIAMETRE,(SIDE_LENGTH/2)-PARTICLE_DIAMETRE)),radius=PARTICLE_DIAMETRE/2,make_trail=MAKE_TRAILS)
        
        #it is assumed that the x,y,z components are NOT normalized to the magnitude.   
        self.speed = _speed
        self.vector = vector(x_component,y_component,z_component).norm() * self.speed
        

    def update_pos(self):
        #change position
        self.pos.x += (1/FRAME_RATE) * self.vector.x
        self.pos.y += (1/FRAME_RATE) * self.vector.y
        self.pos.z += (1/FRAME_RATE) * self.vector.z

        #check for collisions
        if (self.pos.x + self.radius) >= (right_wall.pos.x-(WALL_WIDTH/2)):
            reflect(self,'x')
            self.pos.x = (right_wall.pos.x-(WALL_WIDTH/2))-self.radius
        elif (self.pos.x - self.radius) <= (left_wall.pos.x+(WALL_WIDTH/2)):
            reflect(self,'x')
            self.pos.x = left_wall.pos.x+(WALL_WIDTH/2)+self.radius
        if (self.pos.y + self.radius) >= (top_wall.pos.y-(WALL_WIDTH/2)):
            reflect(self,'y')
            self.pos.y = top_wall.pos.y-(WALL_WIDTH/2)-self.radius
        elif (self.pos.y - self.radius) <= (bottom_wall.pos.y+(WALL_WIDTH/2)):
            reflect(self,'y')
            self.pos.y = bottom_wall.pos.y+(WALL_WIDTH/2)+self.radius
        if (self.pos.z - self.radius) <= (back_wall.pos.z+(WALL_WIDTH/2)):
            reflect(self,'z')
            self.pos.z = back_wall.pos.z+(WALL_WIDTH/2)+self.radius
        elif (self.pos.z + self.radius) >= ((SIDE_LENGTH/2)):
            reflect(self,'z')
            self.pos.z = (SIDE_LENGTH/2)-self.radius
    

NUMBER_OF_TRIALS = 4
MIN_TIME_TO_RECORD = 2 #seconds; time after start of trial where data isn't recorded
TIME = 20 #seconds; total time per trial

#TRIALS
temperatures = [50,100,200] 
num_particles = [4]
volumes = [1000]

if __name__ == '__main__':
    for temp in temperatures:
        for volume_ in volumes:
            for num_particle in num_particles:
                #start new measurement
                TEMPERATURE = temp
                NUMBER_OF_PARTICLES = num_particle
                VOLUME = volume_

                PARTICLE_SPEED = math.pow((3 * BOLTZMANN_CONSTANT * TEMPERATURE)/PARTICLE_MASS,0.5)
                SIDE_LENGTH = math.pow(VOLUME,1/3)#m
                TOTAL_WALL_AREA = 6*math.pow(SIDE_LENGTH,2)#m^2
                
                spawn_walls()
                plt.clf()

                measurement_means = []
                measurement_factors = []
                
                for num_trial in range(NUMBER_OF_TRIALS):
                    #start new trial
                    start_time = time.time()

                    counter = float()
                    global_particle_velocity_change = float()
                    trial_times = list()
                    trial_factors = list()

                    spawn_particles()

                    while counter <= TIME:
                        #start simulation

                        rate(FRAME_RATE)
                        
                        for particle in particles:
                            particle.update_pos()

                            for other_particle in particles:
                                #checking for particle-particle collision
                                if not particle == other_particle:
                                    a = (particle.pos.x,particle.pos.y,particle.pos.z)
                                    b = (other_particle.pos.x,other_particle.pos.y,other_particle.pos.z)
                                    dist = distance.euclidean(a,b)
                                    if dist <= PARTICLE_DIAMETRE:
                                        #two particles have collided!
                                        delta_x = particle.pos.x - other_particle.pos.x
                                        delta_y = particle.pos.y - other_particle.pos.y
                                        delta_z = particle.pos.z - other_particle.pos.z
                                        reflection_vector = vector(delta_x,delta_y,delta_z).norm()
                                        #since collisions are perfectly elastic and their masses are the same, speeds are traded; since velocities are same, magnitude is consistent.
                                        particle.vector = vector(
                                            reflection_vector.x * other_particle.speed, 
                                            reflection_vector.z * other_particle.speed,
                                            reflection_vector.y * other_particle.speed,
                                        )
                                        other_particle.vector = vector(
                                            reflection_vector.x * -particle.speed, 
                                            reflection_vector.z * -particle.speed,
                                            reflection_vector.y * -particle.speed,
                                        )
                                        old_speed = particle.speed
                                        particle.speed = other_particle.speed
                                        other_particle.speed = old_speed
                        
                        #append info to graphs
                        counter = time.time()-start_time
                        if counter > 0:
                            force = (PARTICLE_MASS * global_particle_velocity_change)/counter
                            pressure = force/TOTAL_WALL_AREA
                            expected_pressure = BOLTZMANN_CONSTANT*NUMBER_OF_PARTICLES*TEMPERATURE/VOLUME
                            factor = (pressure*VOLUME)/(NUMBER_OF_PARTICLES*BOLTZMANN_CONSTANT*TEMPERATURE)
                            if counter >= MIN_TIME_TO_RECORD:
                                trial_times.append(counter)
                                trial_factors.append(factor)

                    #finished trial 
                    plt.plot(trial_times,trial_factors)
                    measurement_factors.extend(trial_factors)
                    measurement_means.append(sum(trial_factors)/len(trial_factors))
               
                #finished measurements
                plt.xlabel('Time(seconds)')
                plt.ylabel('Quotient(1 expected)')

                correct_vals = [1 for i in range(len(measurement_factors))]
                rmse = 'RMSE: ' + str(round(mean_squared_error(measurement_factors,correct_vals),4))
                mean = 'Mean: ' + str(round(sum(measurement_means)/NUMBER_OF_TRIALS,4))
                
                graph_title = 'Number of Particles: ' + str(NUMBER_OF_PARTICLES) + ', Temperature: ' + str(TEMPERATURE) + ', Volume: ' + str(VOLUME)
                title = '\n'.join([graph_title,mean,rmse])
                plt.title(title,fontsize=10)
                #this is really lazy code:
                graph_title = graph_title.replace(' ', '')
                graph_title = graph_title.replace(':',"")
                graph_title = graph_title.replace(',','')
                
                
                plt.savefig(graph_title + '.png')
                print('made ' + graph_title + '\n' + rmse)