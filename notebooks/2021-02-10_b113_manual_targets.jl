using LightXML, Random
##
fishDir = "/mnt/deissero/users/tyler/tmp/"
nNeurons = 10
outsideBrainNum = nNeurons+1
nReps = 8
trialOrder = repeat(collect(1:nNeurons),nReps)
trialOrder = trialOrder[randperm(nNeurons*nReps)]

xmlSeries = XMLDocument() # create & attach a root node
xSeriesRoot = create_root(xmlSeries, "PVSavedMarkPointSeriesElements")
set_attribute(xSeriesRoot, "Iterations", "1")
set_attribute(xSeriesRoot, "IterationDelay", "5000.0")


for i in trialOrder

    mpe = new_child(xSeriesRoot, "PVMarkPointElement")
    set_attribute(mpe, "Repetitions", "5")
    set_attribute(mpe, "UncagingLaser", "Uncaging")
    set_attribute(mpe, "UncagingLaserPower", "0")
    set_attribute(mpe, "TriggerFrequency", "None")
    set_attribute(mpe, "TriggerSelection", "None")
    set_attribute(mpe, "TriggerCount", "1")
    set_attribute(mpe, "AsyncSyncFrequency", "None")
    set_attribute(mpe, "VoltageOutputCategoryName", "None")
    set_attribute(mpe, "VoltageRecCategoryName", "Current")
    set_attribute(mpe, "parameterSet", "CurrentSettings")

    # nested under <PVMarkPointElement>!
    gpe = new_child(mpe, "PVGalvoPointElement")
    if i==1
        set_attribute(gpe, "InitialDelay", "15")
    else
        set_attribute(gpe, "InitialDelay", "0.12")
    end
    set_attribute(gpe, "InterPointDelay", "248")
    set_attribute(gpe, "Duration", "170") # from exported file but ????
    # set_attribute(gpe, "Duration", "2")
    set_attribute(gpe, "SpiralRevolutions", "5")
    set_attribute(gpe, "AllPointsAtOnce", "False")
    set_attribute(gpe, "Points", "Point $i")
    set_attribute(gpe, "Indices", "$i" )
end

save_file(xmlSeries, joinpath(fishDir, "2020-02-10_b113_$(nNeurons)cell_sequential.xml"))