function max_95ci(x)
    confint(OneSampleTTest(x))[2]
end

function min_95ci(x)
    confint(OneSampleTTest(x))[1]
end