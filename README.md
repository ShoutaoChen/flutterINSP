# flutterINSP
PROGRAM DESCRIPTION
-------------------
This MATLAB program analyzes the flutter stability of suspension bridges during the erection stage, helping to determine optimal erection start dates and durations while calculating corresponding flutter safety probabilities.

MAIN FEATURES
-------------

Flutter Stability Inspection: Evaluates suspension bridge stability during erection stages

Erection Schedule Optimization: Identifies feasible erection start dates and durations

Probability Analysis: Computes flutter safety probabilities for different erection schemes

Multi-Location Support: Supports analysis for six different cities with varying terrain conditions

CASE BRIDGES
------------

The Severn Bridge in UK

The Xihoumen Bridge in China

CASE LOCATIONS AND TERRAIN
--------------------------
Cities: Dalian, Qingdao, Hangzhou, Taipei, Xiamen, Hong Kong
Terrain Roughness Categories: A, B, C, D

PREREQUISITES
-------------

MATLAB R2022a or compatible version

Historical wind speed data for your bridge location

FILE STRUCTURE
--------------
root/
├── Feasible_start_date.m          # Main program for feasible start date analysis
├── Flutter_safety_probability.m   # Main program for flutter safety probability
├── GAMMA_T.CST                    # Configuration file
└── WDSP_database/                 # Historical wind speed database
    ├── Dalian.txt
    ├── Qingdao.txt
    ├── Hangzhou.txt
    ├── Taipei.txt
    ├── Xiamen.txt
    └── Hong Kong.txt

DATABASE FORMAT
---------------
Each text file contains seven columns:

Column 1: Weather station ID

Columns 2-4: Year, Month, Day

Column 5: Daily maximum 10-minute mean wind speed (knots)

Column 6: Daily maximum 3-second mean wind speed (knots)

Column 7: Daily average wind speed (knots)

USAGE INSTRUCTIONS
------------------

FEASIBLE ERECTION START DATE ANALYSIS
Run Feasible_start_date.m with the total erection duration (d_tot) as input to obtain the feasible timeline of erection start dates for your selected city and terrain category.

FLUTTER SAFETY PROBABILITY ANALYSIS
Run Flutter_safety_probability.m with specified ranges for:

d_tot_range: Range of total erection durations

d_str_range: Range of erection start dates
This generates a cloud map showing flutter safety probabilities for all possible erection schemes.

CUSTOMIZATION FOR YOUR BRIDGE
------------------------------
To apply this program to your specific bridge:

Replace the historical wind speed data with data from your bridge location

Update the flutter critical wind speeds and basic bridge parameters (main span length, deck height) with your specific values

Modify terrain roughness categories as needed for your site

THEORETICAL BACKGROUND
For detailed theoretical foundations, please refer to the paper:
"Feasible Erection Schedules of Suspension Bridge Considering Flutter Safety Under Extreme Wind Speeds"

MATLAB VERSION
--------------
Developed and tested with MATLAB R2022a

CONTACT INFORMATION
-------------------
Email: 244807009@csu.edu.cn

