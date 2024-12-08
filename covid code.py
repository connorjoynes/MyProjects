import random
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import numpy as np
import scipy.stats as stats

#=====================functions===================

def tables (x,y,table_width,table_length,room_size,current_tables):
    '''

    Parameters
    ----------
    x : a random x value within bounds of the pub.
        
    y :  a random y value within bounds of the pub.
    
    table_width : the width of the table.
    
    table_length : length of the table.
        
    room_size : the dimensions of the room.
        
    current_tables : a list of the current tables.
        

    Returns:
        False if the random coordinates overlap with an existing table.
        True if the random coordinated don't overlap any tables.
    '''
    global table_count, social_distancing # allows the function to see these variables
    #checks if social disantcing is turned on
    if social_distancing == True:
        #if it is increase the table distancing to 2m as regulated
        distancing = 2
    else:
        #otherwise keep it as 0.5m
        distancing = 0.5
    #just chacks the table has dimensions
    if table_width > 0 and table_length>0:
        #loops over exsisting tables
        for i in range(len(current_tables)):
            #checks if tables overlap
            if current_tables[i][0]-table_width-distancing < x < current_tables[i][0]+table_width+distancing and current_tables[i][1]-table_length-(distancing/2) < y < current_tables[i][1]+table_length+(distancing/2):
                #if tables do overlap return false
                return False
        #if the new table doesnt overlap with any then begin adding it to the list
        #create a temporary list
        temp = []
        #append the tables x and y coords
        temp.append(x)
        temp.append(y)
        #then add the temp list to the current tables list
        current_tables.append(temp)
        # and increment the number of tables by 1
        table_count = table_count + 1
        return True

def make_table(ax,current_tables, width, height,number_of_tables):
    '''

    Parameters
    ----------
    ax : the current set of axis to add the tables too.
        
    current_tables : list of current tables.
        
    width : table width.
        
    height : table length
        
    number_of_tables : the number of tables
        

    Returns
    -------
    artists : returns the plot of tables.
    

    '''
    
    # Loop over all the tables
    for i in range (number_of_tables):
        #produce a rectangle shape for each table
        table = [Rectangle((current_tables[i][0],current_tables[i][1]), width, height)]

        #this adds the edge of the tables which will stand out more
        pc1 = PatchCollection(table, facecolor='none', alpha=1,
                         edgecolor='darkslategrey',zorder = 3)
        # Create patch collection with specified colour/alpha
        pc = PatchCollection(table, facecolor='darkred', alpha=1,
                          edgecolor='k',zorder = 1)

        # Add collection to axes
        ax.add_collection(pc)
        ax.add_collection(pc1)

        # Plot errorbars
        artists = ax.scatter(current_tables[i][0], current_tables[i][1], s = 0)

    return artists

def customers(current_tables,width,length):
    '''
    
    Parameters
    ----------
    current_tables : list of the tables.
        
    width : width of the table.
        
    length : length of the table.
        

    Returns
    -------
    positions : a list of the peoples positions.       

    '''
    global reduced_mass #allows function to see this variable
    #list will store the positions
    positions = []
    #a number of tables
    vacant_tables = int(len(current_tables))
    #checks if reduced mass is turned on
    if reduced_mass  == True:
        # if it is then make it so half the tables remain empty
        vacant_tables = int(vacant_tables/2)
    #loop over all the tables
    for table in range(vacant_tables):
        #a kind of temporary list for each table
        seats = []
        #loops over all of the avaialble seats
        for i in range(6):
            #left hand side of the table
            if i ==0 or i == 2 or i == 4: 
                x = current_tables[table][0]
                y = current_tables[table][1]+((i/2)+0.5)
                #appends the seat positions to a different temporary list
                temp = []
                temp.append(x)
                temp.append(y)
                #adds the temp list to seats
                seats.append(temp)
            #same but for the right of the table
            if i ==1 or i == 3 or i == 5:
                x = current_tables[table][0]+(width)
                y = current_tables[table][1]+(((i-1)/2)+0.5)
                temp = []
                temp.append(x)
                temp.append(y)
                seats.append(temp)
        #randomly generate how many people are sat at this table        
        sat = random.randint(1,6)
        #loops over how many people are sat
        for j in range(sat):
            #adds each person who is sat from the previous list
            positions.append(seats[j])
    return positions
    

def infection_path(time, spreaders, room_size, current_tables,number_of_people,start,catch_distance):
    '''
    

    Parameters
    ----------
    time : how long the simulation will run for
        
    spreaders : how many initailly contagious people there are
        
    room_size : dimensions of the room
        
    current_tables : list of the tables
        
    number_of_people : number of people in the pub
        
    start : positions of the people in the pub
        
    catch_distance : the range at which people will catch COVID
        

    Returns
    -------
    X : a list of the X coordinates of the virus random walk
        
    Y : a list of the Y coordinates of the virus random walk
        
    patient_zero : list of the people who start ill
        
    infected : list of the pople infected
        
    virus_count : how many times the person will cough and produce a new virus random walk
        
    infection_amount : infection rate of the simulation (%)
        
    vaccinated_person : boolean to say if vaccinations are turned on or not.
        

    '''
    global linger_time, window_y,repeating,info,ventilation#gives the function to these variables
    virus_count = int(time/3600)#coughs once an hour
    
    # X and Y coordinates, is a list for each virus and then a list for each spreader              
    X = [[[] for _ in range(virus_count)] for _ in range(spreaders)]  # List to store X positions for each virus
    Y = [[[] for _ in range(virus_count)] for _ in range(spreaders)]  # List to store Y positions for each virus 
    # List to track infection status of each person
    patient_zero = [False] * number_of_people
    #list of infected people
    person = []
    #keeps a check on who how close people are to being infected.
    infection_check = [0 for _ in range(number_of_people)]
    #variable for loop
    count_spreaders=0
    #loop until all spreaders have been designated a person
    while count_spreaders < spreaders:
        # Randomly select starting people to be infected
        Random = random.randint(0, number_of_people - 1)
        #check the random person isnt already infected
        if Random not in person:
            #add then to the list so we know they are infected
            person.append(Random)
            #assign this person as patient zero
            patient_zero[person[count_spreaders]] = True
            #increment loop counter
            count_spreaders += 1
    #list of vaccinated people       
    vaccinated_person = []
    #if vaccines are turned on
    if vaccinated == True:
        #loops over 93.6% of the population
        for i in range(int(number_of_people*(93.6/100))):
            #makes each person who is looped over vaccinated
            vaccinated_person.append(i)
    #the value that if a random counter reaches the virus path will stop        
    stopping_value = 6
    #checks if ventilation is turned on
    if ventilation == True:
        #if it is set stopping value to 1/3 of the original
        stopping_value = 2
    #set infected to patient zero
    infected = patient_zero
    #keeps count on how many people get newly infected
    infection_count = 0
    #keeps a track on the stopping chance of each virus
    stopped = [0]*virus_count
    #loop over the time
    for current_time in range(time + 1):
        #loop over each virus
        for virus in range(virus_count):
            #checks if the virus has stopped
            if stopped[virus] < stopping_value:  
                #if not then loop over each spreader
                for current_carrier in range(spreaders):
    
                    #if its the first step
                    if current_time < 1:
                        #start at the spreaders location
                        x = start[person[current_carrier]][0]
                        y = start[person[current_carrier]][1]   
                        #set environment to inside the pub
                        outside = False
                        inbounds = True
                    else:
                        #if not the first step then assume virus step is outside bounds
                        inbounds = False
                        #create a random value for cahnce to stop virus
                        stop_chance=random.random()
                        #check if chance is below critical vlaue
                        if stop_chance < 0.000069:
                            #if it is increment stopped for this virus by 1
                            stopped[virus] += 1
                        #begin a loop counter
                        count = 0
                        #loops over new coordinates until a valid one is found
                        while inbounds == False:
                            #increment the counter
                            count += 1
                            #if counter ever reaches 50 then just set new soords as the last coords
                            if count > 50:
                                #set the virus as inside the pub
                                inbounds = True
                                x = X[current_carrier][virus][-1]
                                y = Y[current_carrier][virus][-1]
                                #breaks the loop
                                break
                            #assume the virus stays inside the pub
                            inbounds = True
                            #create a random length to walk
                            length = random.uniform(0, 0.5/(virus+1))
                            #create a random direction to walk
                            direction = 2 * np.pi * random.uniform(0, 1)
                            #add this length and direction to the last set of coords
                            x = length * np.cos(direction) + X[current_carrier][virus][-1]
                            y = length * np.sin(direction) + Y[current_carrier][virus][-1]
     
                        #checks they stay in the room if windows are shut
                            if open_windows == False:
                                if x < 0 or x > room_size:
                                    inbounds = False
                                if y < 0 or y > room_size:
                                    inbounds = False
                        #if windows are open the virus can leave via windows
                            if open_windows == True:
                                #loops over both sets of windows
                                for j in range(int(len(window_y)/2)): 
                                    #checks if the virus is inside the pub
                                    if outside == False:
                                        #checks if the virus is next to the windows
                                        if (X[current_carrier][virus][-1] < 1 or X[current_carrier][virus][-1] > room_size-1) and window_y[0][j] < Y[current_carrier][virus][-1] < window_y[1][j] :                           
                                            #if the new coords are outside the pub then set status to outside
                                            if x < 0 or x > room_size:
                                                outside = True
                                        #checks if virus is not next to the window
                                        else:
                                            #then virus must remain in the room
                                            if x < 0 or x > room_size:
                                                inbounds = False
                                            if y < 0 or y > room_size:
                                                inbounds = False
                                    #checks if the virus is already outside
                                    if outside == True:
                                        #do another random check to see if virus stops
                                        stop_chance=random.random()
                                        #larger chance to imitate wind taking virus away
                                        if stop_chance < 0.0001:
                                            # if succeed then increment stop by 1
                                            stopped[virus] += 1
                                        #checks that the virus is withn the simulation still
                                        if x < -outside_area or x > room_size+outside_area:
                                            inbounds = False
                                        if y < -outside_area or y > room_size+outside_area:
                                            inbounds = False
                                        #checks if the virus is near the window    
                                        if (X[current_carrier][virus][-1] > -1 or X[current_carrier][virus][-1] < room_size+1) and window_y[0][j] < Y[current_carrier][virus][-1] < window_y[1][j] :
                                           #cheks if the virus has come back inside
                                            if 0 < x < room_size:
                                                #set status to inside
                                                outside = False    
                                        else:
                                            #if its not near a window then make sure virus stays outside
                                            if (0<x<room_size and 0<y<room_size):
                                                inbounds = False
  
                    #checks infecting others
                    #loops over every person in the pub
                    for check in range(number_of_people):
                        #checks we arent comparing to the same person and checks if the current person is healthy 
                        if check != person[current_carrier] and infected[check] == False:
                            #checks if the virus is close enough to the person
                            if  start[check][0]-catch_distance < x < start[check][0]+catch_distance and start[check][1]-catch_distance < y < start[check][1]+catch_distance:
                                #checks the person isnt vaccinated
                                if check not in vaccinated_person:
                                    #increment infection status of that person by 1
                                    infection_check[check] += 1 
                                else:
                                    #if vaccinated increment by 0.05 instead
                                    infection_check[check] += 0.05 
                                #checks if the incremented amount of that person has reached the cap
                                if infection_check[check] == infection_number:
                                    #then add them to infected list
                                    infected[check] = True
                                    #increment the number of newly infected by 1
                                    infection_count += 1
                                    #adjust so tracker goes up by 1 from this time on
                                    for second in range(current_time,time):
                                        infection_time[second] = infection_count
                                    break
        
                      
                    #adds the new coords to the list of current coords
                    X[current_carrier][virus].append(x)
                    Y[current_carrier][virus].append(y) 
            # if the virus has stopped
            else:
                #then the new coord is the last current coord
                X[current_carrier][virus].append(X[current_carrier][virus][-1])
                Y[current_carrier][virus].append(Y[current_carrier][virus][-1])
    #checks if reduced mass is turned on
    if reduced_mass == True:
        #assume half the people remain home and healthy
        number_of_people=number_of_people*2
    #calculate the infection rate as a percentage
    infection_amount = (((infection_count+spreaders)/number_of_people)*100)
    #generate information if its the first repetition
    if repeating == False:
        info.append(f'\n\n\n{number_of_people} people in the pub')
        info.append(f'{spreaders} people start with covid')
        info.append(f'\npeople newly infected: {infection_count}')
        info.append(f'percentage leaving with covid: {infection_amount:.3} %')
    return X, Y, patient_zero, infected, virus_count,infection_amount,vaccinated_person


#==================parameters===============================
#not repeating
repeating = False
#set time to 4 hours
time = 14400
#set initially contagious people to 1
spreaders = 1
#this determine how many times the virus must cross a person before they become infected
infection_number = 3
#determine the radius people can catch the virus from
catch_distance = 0.45


#pub
#determine the pubs dimensions
room_size = 20
#table width
table_width = 2
#table length
table_length = 3
#sets how many tables are generated in each pub
number_of_tables = 8
#keeps track of how many tables have been generated
table_count = 0
#list of the generated tables
current_tables = []
#sets the range of the outside area
outside_area = 3
#sets the coordinates of the windows
window_x=[[0,0,room_size,room_size],[0,0,room_size,room_size]]
window_y=[[(room_size/5),(3*room_size/5),(room_size/5),(3*room_size/5)],[((room_size/5)+room_size/10),((3*room_size/5)+room_size/10),((room_size/5)+room_size/10),((3*room_size/5)+room_size/10)]]
#list that information is appened to so its seen at the end
info = []

#precautionary measures
social_distancing = True
social_distancing = False
reduced_mass = True
reduced_mass = False
masks = True
masks = False
open_windows = True
open_windows = False
ventilation = True
ventilation = False
reduced_time = True
reduced_time = False
vaccinated = True
vaccinated = False
#checks if reduced time is turned on
if reduced_time == True:
    #just reduces the time to 2 hours
    time = int(time/2)
#keeps a track at when people get infected
infection_time = [0 for i in range(time)]

#checks if masks are turned on
if masks == True:
    #reduces infection by 90% by increasing the required number of times with covid to increase by 10
    infection_number = infection_number*10
    # reduces the catch distance
    catch_distance = 0.18
    
#=========================code==================================
# sets a loop counter
count = 0
#loops until all tables have been placed
while table_count != number_of_tables:
    #increments the counter by 1
    count += 1
    #runs the table function
    tables(random.randrange(1,room_size-table_width,1),random.randrange(1,room_size-table_length,1),table_width,table_length,room_size,current_tables)
    #if cant place all tables after 100 loops then reduce counter by 1 and try again
    if count > 100:
        table_count -= 1

#runs the function that assigns poeple to the pub
people = customers(current_tables, table_width, table_length)
#variable that records how mnay people there are
number_of_people = len(people)
#runs the function that does the random walk
X, Y, patient_zero, infected, virus_count, infection_count,vaccinated_people = infection_path(time, spreaders, room_size, current_tables,number_of_people,people,catch_distance)

#==================simulation=====================
#set 2 subplots side by side
fig,(ax1,ax2)=plt.subplots(1,2,figsize = (15,10))
#generates the tables
_ = make_table(ax1,current_tables, table_width, table_length,number_of_tables)
#seperates the list of the poples coords into x and y values
sitting_x = [x for x,y in people]
sitting_y = [y for x,y in people]
#loops over all people
for i in range(number_of_people):
    #plots all vaccinated pople ad cyan
    if i in vaccinated_people:
        ax1.scatter(sitting_x[i],sitting_y[i],zorder = 4, c = 'cyan')
    #plots all newly infected people as orange
    elif infected[i]:
        ax1.scatter(sitting_x[i],sitting_y[i],zorder = 4, c = 'orange')
    # plots all the healthy people as green
    else:
        ax1.scatter(sitting_x[i],sitting_y[i],zorder = 4, c = 'limegreen')

#loops over all spreaders
for i in range(spreaders):
    for j in range(virus_count):
        #makes the contagious person indigo
        ax1.scatter(X[i][j][0], Y[i][j][0], color='Indigo', zorder=5)
        # plots the virus random walk
        ax1.scatter(X[i][j], Y[i][j], label=f"{spreaders}",alpha=0.1,s=5,zorder=2,c='red')
#set out the plot with labels and dimensions    
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_title(f'pub after {time/60} minutes')
ax1.set_xlim(-outside_area,room_size+outside_area)
ax1.set_ylim(-outside_area,room_size+outside_area)
#visuals for the pub
#if turned open windows on then make windows more obvious
if open_windows == True:
    window_visual = 5
else:
    window_visual = 3
    #create a door
door_x = 4
door_y = -0.6
#plots all the windows
ax1.plot(window_x,window_y,c = 'azure',linewidth=window_visual,zorder = 4)
#plots the door
ax1.plot([8,door_x],[0,door_y],c='brown',linewidth = 3,zorder = 4)
#plot the grass
ax1.fill_between([-outside_area,room_size+outside_area],-outside_area,room_size+outside_area,zorder=0,color = 'forestgreen', alpha=0.75)
#plots the pub floor
ax1.fill_between([0,room_size],0,room_size,zorder=0,color = 'peru', alpha=0.75)
#define a percentage size of the pub 
length = outside_area/(room_size+outside_area)
#create darker boarder outlines of the pub walls
ax1.axhline(0,length,(1-length),c='darkslategrey',zorder=3)
ax1.axhline(room_size,length,(1-length),c='darkslategrey',zorder=3)
ax1.axvline(0,length,(1-length),c='darkslategrey',zorder=3)
ax1.axvline(room_size,length,(1-length),c='darkslategrey',zorder=3)

#==============graph=============================
#create an array for time which can be plotted
timing = np.linspace(0,time/60,time)
#plot how the number of infected increases over time
plt.plot(timing,infection_time)
#set graph up for this plot
ax2.set_xlabel('time (min)')
ax2.set_ylabel('newly infected', labelpad = 0)
ax2.set_title('infection per minute')
#shows the top bound
ax2.axhline(number_of_people-spreaders,linestyle = '--', color = 'grey')
ax2.grid()
ax2.set_xlim(0,time/60)
ax2.set_ylim(0,number_of_people+1-spreaders)
plt.show()

# ===================repeat simulations==========================
#now we are repeating it all again
repeating = True
#sets repetitions too 100
repetitions = 100
#keeps track of each repetitions infection rate
percentage_infected = []
#loops over for every repetions
for repeat in range(repetitions):
    #shows how far along we are
    print(f'percentage complete {repeat}%')
    #sets a loop counter
    count = 0
    #list of all current tables
    current_tables = []
    #keeps count of the number of current tables
    table_count = 0
    #loops until all tables are generated
    while table_count != number_of_tables:
        #increment loop counter by 1
        count += 1
        #run table function to generate tables
        tables(random.randrange(1,room_size-table_width,1),random.randrange(1,room_size-table_length,1),table_width,table_length,room_size,current_tables)
        #same as before this just stops infinite loops
        if count > 100:
            table_count -= 1
    #generates the people in the pub        
    people = customers(current_tables, table_width, table_length)
    #how many people are in the pub
    number_of_people = len(people)
    #runs the random walk function
    X, Y, patient_zero, infected, virus_count, infection_count,vaccinated_people = infection_path(time, spreaders, room_size, current_tables,number_of_people,people,catch_distance)
    #adds the infection rate of this simulation to the list    
    percentage_infected.append(infection_count)
#finds the average of all the simulation's infections rates
infection_average=sum(percentage_infected)/len(percentage_infected)

#=====plot averages===========
#creates an x axis of repetitions
r = np.linspace(1,repetitions,repetitions)
#plots the simulation's infection rates
plt.plot(r,percentage_infected, label = 'infected')
#adds the average to this plot
plt.axhline(infection_average,c='orange',label = f'average infection ({infection_average:.3}%)')
#sets the graph for the plot
plt.legend()
plt.grid()
plt.ylim(0,101)
plt.xlabel('number of repetitions')
plt.ylabel('percentage infected (%)')
plt.show()
#loops over all the information generated
for _ in range(len(info)):
    # and prints it all
    print(info[_])
#also prints the average infected rate over the 100 repetitions
print('average percentage infected over 100 repetitions:',infection_average)

#======================erros============================
#calculates the standard deviation
std = np.std(percentage_infected)
#prints it
print('standard deviation:',std)

#creates an array of x values for gaussian plot
x = np.linspace(infection_average - 3*std,100, 100)
#creates the gaussian plot
y=stats.norm.pdf(x, infection_average, std)
#plots the gaussian fit
plt.plot(x, y,color='orange')
#plots the average
plt.axvline(infection_average,c='red', linestyle = '--', label = f'Average {infection_average:.3}%')
#set up for histogram
num_bins = 7
bin_edges = np.linspace(50, 103, num_bins + 1)
binned_values = np.histogram(percentage_infected, bins=bin_edges)[0]
#plots histogram
plt.bar(bin_edges[:-1], binned_values/858, width=np.diff(bin_edges), edgecolor='black')
plt.xlabel('Percentage Infected (%)')
plt.ylabel('Probability Density')
plt.legend()
plt.show()
#calculate and print the standard error
SE = std / np.sqrt(repetitions)
print('standard error:',SE)

#=====================other plots==========================
#a list of all our precautions
name = [f'no\n precautions',f'social\n distancing',f'reduced\n mass','masks',f'open \nwindows','ventilation',f'reduced \ntime', 'vaccinated']
#a list of all their averages over 100 cycles
av = [84.2431,83.543, 41.86,25.87,68.7902,80.7266,61.631,9.9066]
#sets figures size
plt.figure(figsize=(15,10))
#create the bar chart of the individual precautions infections rates
plt.bar(name,av,color='forestgreen')
plt.xticks(size=17)
plt.grid()
plt.ylim(0,100)
plt.ylabel('Infected (%)',fontsize = 'xx-large')
plt.show()

#this will take each set of precautions average and STD and plot the gaussian fits
x = np.linspace(0, 100, 100)
infection_average = 84.2
std=11.1
y=stats.norm.pdf(x, infection_average, std)
plt.plot(x, y,label='no precautions (84.2%)')
infection_average = 44.5
std=6.4
y=stats.norm.pdf(x, infection_average, std)
plt.plot(x,y,label='precaution set 1 (44.5%)')
infection_average = 15.1
std=6.2
y=stats.norm.pdf(x, infection_average, std)
plt.plot(x, y,label='precaution set 2 (15.1%)')
infection_average = 5.51
std=2.8
y=stats.norm.pdf(x, infection_average, std)
plt.plot(x, y,label='precaution set 3 (5.51%)')

plt.axvline(100,c='dimgrey')
plt.axvline(0,c='dimgrey')
plt.grid()
plt.legend(loc = 'upper right')
plt.ylabel('Probability Density')
plt.xlabel('Infection Rate (%)')
plt.show()
