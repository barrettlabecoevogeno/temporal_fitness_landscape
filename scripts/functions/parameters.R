# # # # # # # # # # # # # # # # # # # # 
# Author: Marc-Olivier Beausoleil 
# Date: 12 March 2020
# McGill University 
# List of parameters to load
# # # # # # # # # # # # # # # # # # # # 

# Select the variables to run the scripts ---------------------------------
save.data = "./output/biotic.factors.on.survival_Andrew_meeting_changing_PCA_SCORE_for_only_fortis.RData"

pdf.output <- FALSE
find.peaks.and.valleys = FALSE
standard.ized = TRUE
# If you want only the linear βX not the βx and γX^2
linear.only = TRUE  
# If you want orthogonal X values: BEWARE! This is not the same as modeling an ORTHOGONAL regression 
orthogonal.x = FALSE

# prepare the variables that I'm going to record while iterating 
fit.grad.table = list()
model.list.logistic = NULL
pseudr = NULL
gofitfit = NULL
glm.model.list = NULL
yearly.number.of.id = NULL
all.ranges.x = NULL
list.min.fit = NULL
list.min.fit.trait = NULL
original.x = NULL
midd.list = NULL
my.eco.evo.df = NULL
full.data = NULL
final.df.new.analysis = NULL
newlist = list()
oldxlist = list()
oldzlist = list()
old_beak_L_list = list()
old_beak_W_list = list()
old_beak_D_list = list()
old_band_list = list()
model.list = list()
survived.list = list()
raw.data = NULL
model.list1 = NULL
model.list2 = NULL
