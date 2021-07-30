using FileIO, Unitful, AxisArrays, StaticArrays,
    LibExpat, Glob
import Base.Threads.@threads
using Unitful: μm, s

"Return units in microns from Bruker xml file."
function read_microns_per_pixel(xml_file::String)
    read_microns_per_pixel(read_xml(xml_file))
end

function read_microns_per_pixel(xml::ETree)    
    micronsPerPixel_xml = xml[xpath"""//PVStateValue[@key="micronsPerPixel"]"""][1]
    # parse xml
    lookup_μm = axis -> etree -> parse(Float64,
        micronsPerPixel_xml[xpath"""IndexedValue[@index="$axis"]"""][1].attr["value"])μm
    microscope_units = funprod(map(lookup_μm, ["YAxis", "XAxis", "ZAxis"])...)(micronsPerPixel_xml)
    # zseries = centered(zseries) # need centered for qd registration
    microscope_units
end

"Return z-plane of Z or ETL for first DataFrame in sequence."
function read_first_zaxis(xml::ETree, device="Z")
    parse(Float64, xml[xpath"""/PVScan/Sequence/Frame[1]/PVStateShard/PVStateValue[@key="positionCurrent"]/SubindexedValues[@index="ZAxis"]/SubindexedValue[@description="$device"]/@value"""][1])
end

"Return z-plane of Z or ETL for all frames in sequence."
function read_all_zaxis(xml::ETree, device="Z")
    parse.(Float64, xml[xpath"""/PVScan/Sequence/Frame/PVStateShard/PVStateValue[@key="positionCurrent"]/SubindexedValues[@index="ZAxis"]/SubindexedValue[@description="$device"]/@value"""])
end

