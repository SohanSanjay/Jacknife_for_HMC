# Import header files
import numpy as np
import csv

#Parameters of the lattice simulation

L = 8; # size of the lattice : will be usefull when dealing with bond susceptability
       # and also to find the total number of values in each file (for charge susceptability)

N_bins = 100; # Number of bins that are used 

#Each simulation result were written in to several different files
#Each file were save in the same local folder with number of itteration added to the back of the file name.

File_Path = "C:/Users/sohan/Downloads/beta_18/ssh_hmc_square-1/PairSusc_momentum_f/"; 
File_Name = "PairSusc_momentum_";

#For each differnt file read the seond line and the second column (q=(0,0))

#index : the number that added to the end of the file name
#Index of file : index number of the first coloumn in the data output files
#File_name : name of the file
#File_path : path of the file
#read_col : the coloumn of the data value thats going to read

#for pair_susceptability index_of_file = 1 and read_col = 1
#for bond_pair susceptability need to run the csv_read four time and take the sum of the values
#for charge susceptibility nedd to run many time as number of values in the file to take maximum



################################################################################################################################
# The below function is used for all functions that used to read data
# mainly for finding pair susceptability, bond pair susceptability and charge susceptability related readings
################################################################################################################################

def csv_read(index,File_name,File_path,index_of_file,read_col):
    #Full file name and path with file name of the file that need to be read
    File_name_full = File_name+str(index).rjust(5,'0')+str(".out");
    File_read_full = File_path+File_name_full;
    
    #Values that are need to retun
    
    #Read the csv file
    with open(File_read_full, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            row_num = str(row).strip(']').strip('[').split(' ')[0].strip(' \' ');
            χ_val = str(row).strip(']').strip('[').split(' ')[read_col].strip(' \' ');
            if row_num == str(index_of_file):
                return np.float64(χ_val);
                break;
                
################################################################################################################################
# The below function is used to read Pair_susceptability values from raw data files
# index is the added number at the end of the each data file
################################################################################################################################
                
def Pair_susceptibility_data_files(index,File_path):
    File_Name = "PairSusc_momentum_";# name of the file
    index_of_file = 1; #q = (0,0) points 
    read_col = 1;
    χ = csv_read(index,File_Name,File_path,index_of_file,read_col);
    return χ;

################################################################################################################################
# Right another similar function for bond susceptability
# Need to add all four different values to get total bond pair susceptability 
################################################################################################################################




################################################################################################################################
# Right another function for charge_susceptability
# Unlike previous case need to know how many indexes were there in each file
# First read the stat_out file that proces by the previous simulation code and find the index of the maximum obtained value their
# Read each file for similar index and then apply the jacknife method their
# Mean should be similr but the jacknife error should be different
################################################################################################################################






################################################################################################################################
# The below function will calculate mean and jacknife uncertinity for each variable
# This will use the 'Pair_susceptibility_data_files' function
#
# 1. calculate N_bins amout of mean values
# 2. 
################################################################################################################################

def Pair_Jacknife(File_Path):
    x_values = [] #save all N_bins ammout of data values that belongs to the pair_sus.
    x_ex_means = [] #save N_bins ammount of mean values after excluding the one values from the set
    
    # Loop that used to store the pair susceptablity data in x_values array 
    for index in range(1,N_bins+1):
        X_val = Pair_susceptibility_data_files(index,File_Path);
        x_values.append(X_val);
    
    #calculate the N_bins amount of average values and store them on the x_ex_means array
    #First calculate the total sum and exclude each value from the data set and obtain N_bins ammount of average values
    x_total_sum = np.sum(x_values);

    for index in range(1,N_bins+1):
        x_sub_total = x_total_sum - x_values[index-1]; # here use index-1 due to array been on 0 to 99 range
        x_sub_mean = x_sub_total/(N_bins-1);
        x_ex_means.append(x_sub_mean);
    
    #Average value obtained by applying Jacknife method 
    X_avg_JN = np.sum(x_ex_means)/N_bins;

    #Standard diviation obtained by using Jacknife method
    #First subtract the Jacknife mean value from all the values in the array
    x_sub_mean_square = [(x - X_avg_JN)*(x - X_avg_JN) for x in x_ex_means];
    #Second calculate the Standard diviation from this dataset
    X_std_deviation = np.sqrt((N_bins-1)*(np.sum(x_sub_mean_square))/(N_bins));

################################################################################
# For third part of the 
    return X_avg_JN,X_std_deviation;

################################################################################
# The following part contain modifications done to the jacknife method
# 1. First the value is calculated according to the usual jackknife procedure
# 2. The mean and the uncertinity of the mean was calculated using the function before
# 3. Then the modified function for jackknife is run
#       #1 for this function some values were removed fromm the original data set
#       #2 the values removed depending on the mean and the uncertinity
#       #3 the 'in_cons' value is used as an input for the modified function
#       #4 (in_cons*uncertinity) values away from the mean was removed and new data set were created
#       #5 Jackknife method was applied over the new data set to obtain new man and uncertinity value
################################################################################

def Modified_data_JN_Pair(in_cons,File_Path):
    val_avg,val_std = Pair_Jacknife(File_Path);

    #######################################################################################
    #The following content was coppied over frm=om the jacknife function written above
    x_values = [] #save all N_bins ammout of data values that belongs to the pair_sus.
    x_ex_means = [] #save N_bins ammount of mean values after excluding the one values from the set
    
    # Loop that used to store the pair susceptablity data in x_values array 
    for index in range(1,N_bins+1):
        X_val = Pair_susceptibility_data_files(index,File_Path);
        x_values.append(X_val);
    #######################################################################################
    # 'val_std*np.sqrt(N_bins)' is taken insted of 'val_std' because the previous value represents 
    # The uncertinity of the mean

    val_max_tolarance = val_avg+val_std*np.sqrt(N_bins)*in_cons; # values greater than this in 'x_values' are removed
    val_min_tolarance = val_avg-val_std*np.sqrt(N_bins)*in_cons; # values smaller than this in 'x_values' are removed 

    x_modified_values = [];

    for index in range(1,N_bins+1):
        X_val = x_values[index-1];
        if (X_val > val_min_tolarance) & (X_val < val_max_tolarance):
            x_modified_values.append(X_val);
    
    #######################################################################################
    #   #1 Find the number of elements in the 'x_modified_values'  
    #   #2 Apply the jackknife method again over 'x_modified_values'

    N_bins_new = np.array(x_modified_values).size;

    # Following part of the code is a copy of the previous jackknife

    #calculate the N_bins amount of average values and store them on the x_ex_means array
    #First calculate the total sum and exclude each value from the data set and obtain N_bins ammount of average values
    x_total_sum = np.sum(x_modified_values);

    for index in range(1,N_bins_new+1):
        x_sub_total = x_total_sum - x_modified_values[index-1]; # here use index-1 due to array been on 0 to 99 range
        x_sub_mean = x_sub_total/(N_bins_new-1);
        x_ex_means.append(x_sub_mean);
    
    #Average value obtained by applying Jacknife method 
    X_avg_JN = np.sum(x_ex_means)/N_bins_new;

    #Standard diviation obtained by using Jacknife method
    #First subtract the Jacknife mean value from all the values in the array
    x_sub_mean_square = [(x - X_avg_JN)*(x - X_avg_JN) for x in x_ex_means];
    #Second calculate the Standard diviation from this dataset
    X_std_deviation = np.sqrt((N_bins_new-1)*(np.sum(x_sub_mean_square))/(N_bins_new));
    print(x_modified_values)
################################################################################
# For third part of the 
    return X_avg_JN,X_std_deviation;




##############################################################
#This section of the code is used to check the function

val_avg,val_std = Pair_Jacknife(File_Path);
val_avg1,val_std1 = Modified_data_JN_Pair(5,File_Path)
print(val_avg);
print(val_std);
print(val_avg1);
print(val_std1);


