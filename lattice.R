######################################################
### Numerical Simulation of the Ising model        ###
### by Markus Kirkines                             ###
### WT 2017 University of Bonn                     ###
### Computational Physics Course                   ###
######################################################


### set numerical and physical parameters

mc_steps <- 1500   ### 1 MC-Step is N^2 (possible spin-flip) steps
dim_N <- 50        ### dimension of the lattice N^2
J <- 1             ### coupling strength of Hamiltonian (J>0: ferromagnetic behaviour)
k_b <- 1           ### boltzmann constant
T_start <- 1.4     ### start temperature
T_end <- 3.2       ### end temperature
step_size <- 0.025 ### temperature step
iterations <- (T_end-T_start)/step_size ### number of iterations of the main loop
beta <- 1/(k_b*T_start) ### reciprocal temperature

### This function calculates the value of the Energy
### for a given spin.
### The boundary conditions for the nearest neighbours 
### are implemented via the differt if-conditions.
get_spin_energy <- function(mat,i,j){ 
	if(!is.matrix(mat)){
		stop("Argument must be a matrix")
	} else if (i==1 && j==1){
		x=(mat[i+1,j]+mat[i,j+1]+mat[dim_N,j]+mat[i,dim_N])
	} else if (i==dim_N && j==dim_N){
		x=(mat[i-1,j]+mat[i,j-1]+mat[1,j]+mat[i,1])
	} else if (i==dim_N && j==1){
		x=(mat[i-1,j]+mat[i,j+1]+mat[1,j]+mat[i,dim_N])
	} else if (i==1 && j==dim_N){
		x=(mat[i+1,j]+mat[i,j-1]+mat[i,1]+mat[dim_N,j])
	} else if (i==1){
		x=(mat[i+1,j]+mat[i,j+1]+mat[i,j-1]+mat[dim_N,j])
	} else if (i==dim_N){
		x=(mat[i-1,j]+mat[i,j+1]+mat[i,j-1]+mat[1,j])
	} else if (j==1){
		x=(mat[i+1,j]+mat[i-1,j]+mat[i,j+1]+mat[i,dim_N])
	} else if (j==dim_N){
		x=(mat[i+1,j]+mat[i-1,j]+mat[i,j-1]+mat[i,1])
	} else {
		x=(mat[i+1,j]+mat[i-1,j]+mat[i,j-1]+mat[i,j+1])
	}
	return(J*mat[i,j]*x)
}

### This function calculates the total energy 
### of the spin configuration.
### All spin energies are added up and then divided by 4 
### (because of 4 nearest neighbours) and by the number of
### spin particles N².
get_total_energy <-function(mat){
  total <- 0
  for(i in 1:dim_N){
    for(j in 1:dim_N){
      total = total - get_spin_energy(mat,i,j)
    }
  }
  return(total/dim_N^2/4)
}

### This function calculates the absolute value of the 
### magnetisation per spin.
### All the spins are added up and the sum is divided
### by the total number of spin particles N².
get_mag_per_s <- function(mat){
  mag <- 0
  mag = sum(mat)
  mag = abs(mag)
  return(1/dim_N^2*mag)
}

### Hot Spin configuration:
### This function randomly sets up the spins in the lattice with
### values either 1 or -1
init_spin_config_hot <- function(){
  matrix(sample(c(-1, 1), dim_N^2, replace = TRUE), ncol = dim_N, nrow = dim_N)
}

### Cold Spin configuration:
### This function sets up the spins in the lattice all 
### facing in one direction.
init_spin_config_cold <- function(){
  matrix(-1,ncol = dim_N, nrow = dim_N)
}

####################
### MAIN PROGRAM ###
####################
lattice <- array(rep(0,dim_N^2), dim=c(dim_N,dim_N)) ### initialize matrix
lattice <- init_spin_config_cold() ### write all alligned spins into matrix
mag_list <- rep(0,iterations) ### create lists to safe magnetization
ene_list <- rep(0,iterations) ### and energies.
x <- (1:nrow(lattice)) ### put the lattice dimensions as an array
y <- (1:ncol(lattice)) ### for better plotting

### Print the image of the initialized lattice
image(x,y,lattice)
ptm <- proc.time() ### The current time is saved in ptm
pb <- txtProgressBar(min = 0, max = iterations, style = 3) ### set up a progress bar
for (tmp in 1:iterations){ ### MAIN LOOP (iterating over different temperatures)
  T_current <- T_start+(tmp-1)*step_size ### calculate the current temperature
  beta <- 1/(k_b*T_current) ### update beta
  for (step in 1:mc_steps){ ### simulation loop for a single temperature
    for (count in 1:dim_N^2){ ### each step runs N² times
      i <- round(sample(1:dim_N,1)) ### get random coordinates
      j <- round(sample(1:dim_N,1))
      delta_E<-2*get_spin_energy(lattice,i,j) ### calculate possible energy difference
      if (delta_E < 0){ ### check if spin has to be flipped
        lattice[i,j]=-lattice[i,j]
      } else if (exp(-delta_E*beta) >= runif(n=1,min=0,max=1)) {
        lattice[i,j]=-lattice[i,j]
      }
    }
    if (step %% 50 == 0){
      image(x,y,lattice) ### show current lattice
    }
  }
  mag_list[tmp]=get_mag_per_s(lattice) ### save current maginization
  ene_list[tmp]=get_total_energy(lattice) ### and energy to their lists
  lattice <- init_spin_config_cold(lattice) ### reinitialize lattice for the next run
  setTxtProgressBar(pb, tmp) ### update progressbar
}
close(pb) ### close progressbar
print(proc.time() - ptm) ### get elapsed time and print it

### Plot Magnetization and Energy vs Temperaure and write data to file
ttemp <- (1:iterations)
ttemp = T_start+(ttemp-1)*step_size
plot(ttemp,mag_list, xlab = "Temperature",ylab = "Magnetization |m|")
plot(ttemp,ene_list, xlab = "Temperature",ylab = "Total Energy")
df <- data.frame(ttemp,mag_list,ene_list)
write.table(df, "mydata.txt", sep="\t") 