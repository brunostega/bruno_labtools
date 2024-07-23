#!/bin/bash
#SBATCH --job-name=anp_run_300K
#SBATCH -A IscrC_PLIMGO
#SBATCH -p boost_usr_prod
#SBATCH --qos=boost_qos_dbg # uncomment to run in  debug mode
#SBATCH --time 24:00:00     # format: HH:MM:SS, max 24 hours
#SBATCH -N 1                # 1 node
#SBATCH --ntasks-per-node=1 # 4 mpi task per node (to use all four GPUs)
#SBATCH --cpus-per-task=32  # 8 openmp threads per mpi tast (total 32)
#SBATCH --gres=gpu:1        # 4 gpus per node out of 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bruno.stegani@studenti.unimi.it
#SBATCH --exclusive


# Set some environment variables 
FREE_ENERGY=`pwd`
echo "Free energy home directory set to $FREE_ENERGY"
MDP=$FREE_ENERGY/mdps
echo ".mdp files are stored in $MDP"

# Change to the location of your GROMACS-2018 installation
NTOMP=6
PIN=0

for (( i=0; i<26; i++ ))
do

    LAMBDA=$i

    # A new directory will be created for each value of lambda and
    # at each step in the workflow for maximum organization.
    mkdir Lambda_$LAMBDA
    cd Lambda_$LAMBDA

    ##############################
    # ENERGY MINIMIZATION STEEP  #
    ##############################
    echo "Starting minimization for lambda = $LAMBDA..." 

    mkdir EM
    cd EM

    # Iterative calls to grompp and mdrun to run the simulations
    if test -f "min$LAMBDA.gro";then
	   echo "already finished"
    else
 	if [ -f min$LAMBDA.cpt ]; then	    
    		gmx_mpi mdrun -deffnm min$LAMBDA -cpi min$LAMBDA -ntomp $NTOMP -pin on 
	else
		gmx_mpi grompp -f $MDP/EM/ff_em_$LAMBDA.mdp -c $FREE_ENERGY/box.gro -p $FREE_ENERGY/topol_GRETA.top -o min$LAMBDA

    		gmx_mpi mdrun -deffnm min$LAMBDA  -ntomp $NTOMP -pin on -pinoffset $PIN 
    	fi 
    fi
    #####################
    # NVT EQUILIBRATION #
    #####################
    echo "Starting constant volume equilibration..."

    cd ../
    mkdir EQ
    cd EQ

    if test -f "eq$LAMBDA.gro";then
	   echo "already finished"
    else	   
 	if [ -f eq$LAMBDA.cpt ]; then
    		gmx_mpi mdrun -deffnm eq$LAMBDA -cpi eq$LAMBDA -ntomp $NTOMP -pin on  
	else
		gmx_mpi grompp -f $MDP/EQ/ff_aa_eq_$LAMBDA.mdp -c ../EM/min$LAMBDA.gro -p $FREE_ENERGY/topol_GRETA.top -o eq$LAMBDA
	
    		gmx_mpi mdrun -deffnm eq$LAMBDA -ntomp $NTOMP -pin on  -pinoffset $PIN
    	fi 
    	echo "Constant volume equilibration complete."
    	sleep 10
    fi

    #################
    # PRODUCTION MD #
    #################
    echo "Starting production MD simulation..."

    cd ../
    mkdir MD
    cd MD
    if test -f "md$LAMBDA.gro";then
	   echo "already finished"
    else	   
 	if [ -f md$LAMBDA.cpt ];then
    		gmx_mpi mdrun -deffnm md$LAMBDA -cpi md$LAMBDA -ntomp $NTOMP -pin on 
	else 
		gmx_mpi grompp -f $MDP/MD/ff_aa_$LAMBDA.mdp -c ../EQ/eq$LAMBDA.gro -p $FREE_ENERGY/topol_GRETA.top  -o md$LAMBDA.tpr 

    		gmx_mpi mdrun -deffnm md$LAMBDA -ntomp $NTOMP -pin on -pinoffset $PIN
	fi
    	echo "Production MD complete."
    	# End
    	echo "Ending. Job completed for lambda = $LAMBDA"
    fi 
    cd $FREE_ENERGY
done

exit;
