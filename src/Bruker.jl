using FileIO, Unitful, AxisArrays, StaticArrays,
    LibExpat, Glob
import Base.Threads.@threads
using Unitful: μm, s

"Return units in microns from Bruker xml file."
function read_microns_per_pixel(xml_file)
    open(glob("*.xml",xml_file)[1], "r") do io
        xml = read(io, String)
        global xml = xp_parse(xml)
    end;
    
    "apply each function to args"
    funprod(functions...) = (args...)->map(x->x(args...), functions)
    zeroToOne(x) = (x.-minimum(x))/(minimum(x) + maximum(x))
    
    micronsPerPixel_xml = xml[xpath"""//PVStateValue[@key="micronsPerPixel"]"""][1]
    # parse xml
    lookup_μm = axis -> etree -> parse(Float64,
        micronsPerPixel_xml[xpath"""IndexedValue[@index="$axis"]"""][1].attr["value"])μm
    microscope_units = funprod(map(lookup_μm, ["YAxis", "XAxis", "ZAxis"])...)(micronsPerPixel_xml)
    # zseries = centered(zseries) # need centered for qd registration
    microscope_units
end