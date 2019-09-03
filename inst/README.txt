[1]  Call_Model
	[A] Revised_Compile.R - 
		Compiles the individual results files from fitting the model
                in subsets of the full data. (Real Data)

	[B] Revised_Sim_Mix_Calls_MI.R - 
		Killdevil (LSF) submission functions for model fitting in the mixture samples.

	[C] Revised_Sim_MixCode_SI.R - 
		Processes input from Mix_Calls_MI.R and fits the IsoDeconv_MM Model.

	[D] Revised_Sim_PS_calls_new.R
		Killdevil (LSF) submission functions for model fitting in the mixture samples.

	[E] Revised_Sim_PS_codeNew.R -
		Processes input from PS_calls_new.R and fits the IsoDeconv_MM model for pure samples.

[2]  Cuffidff_Out:
	Output from Cuffdiff calls used to generate the Half-Simulated Data.

[3]  Mixture_Creation:
	Programs used to generate the half-simulated mixture data.

	[A] downsample_prog.R - 
		Downsamples the real datasets in order to produce files like (cd4_20 and cd8_80)
		These files will then be pushed together to produce the full mixtures.

	[B] EffectiveLength_Total.R -
		Borrows the functions from Dr. Sun's IsoDot package for computing effective length.

	[C] GeneModel_Build.R -
		Calls the geneModel_multcell_edit.R function to build up the list structure utilized
		by Dr. Sun's IsoDot package.

	[D] geneModel_multcell_edit.R
		Edited versions of Dr. Sun's genemodel() function from the IsoDot package. This function
		is edited to allow for multiple cell types within a single list structure.

	[E] EnsemblIds2Use.RData:
		
	[F] loadData_edit.R
		Support function borrowed from Dr. Sun's IsoDot package to appropriately read in
		data and construct the input list.

	[G] merge_downsamplev2.R
		merges the downsampled files from [A] into composite mixtures. 

	[H] pdist_gen.R:
		Support function borrowed from Dr. Sun's IsoDot package to construct the effective
		length distribution.

	[I] RecharacterizeOutput.R
		Function to be run after geneModel_build.R [C] to recharacterize the output in a manner
		usable by IsoDeconv_MM.

	[J] rem_clust.R
		Preprocessing function designed to remove lowly expressing genes.

	[K] Restricting_Gene_ids_Pairwise.R
		Processes the cuffdiff output to get the analysis genes to use.

[4]  Original Programs
	An Older batch of simulation functions. May not be of use, but included just in case.
	These functions were developed before the multiple initial points (MI) feature was 
	included.

[5]  Process Real Data:
	Files and programs used to process the real data for half simulated mixtures. Download
	.sras, map with star, prep for cufftools and process/clean files.

[6]  Results
	Results of various simulated fits. PlotSpecs_New.R contains the code for visualizing the
	output. Simulation_Function.R contains some support functions for creating the fully
	simulated data. sub_gopt.RData contains a real data derived gene model structure that is
	used for simulation. This ensures reasonable accuracy in the utilized gene models when
	conducting simulations.

[7]  Production Functions (Err_Control_MI_v2).R
	This houses the IsoDeconv_MM Fit functions for mixture data. The code is well documented
	with full functional descriptions at the bottom of the program.

[8]  PS Production Functions (Err Control).R
	This houses the IsoDeconv_MM Fit functions for pure data. The functions in this code are relatively
	simple to understand. Structure of functions and input are very similar to [7].

QUESTIONS: Please email douglas.roy.wilson@gmail.com