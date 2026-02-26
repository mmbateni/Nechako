{standard Water Survey of Canada (WSC) data quality code} While specific definitions can vary slightly by station processing, the standard HYDAT meanings are:
NA (in the Symbol column): In the context of the R dataframe output, this usually means no flag was applied. It indicates standard, approved data with no special conditions (equivalent to a blank space in raw WSC tables).
B: Ice affected. The water level or flow was influenced by ice cover, which may reduce the accuracy of the rating curve.
D: Discharge measured. A field measurement of discharge was taken on this day. This often indicates a shift in the rating curve was applied based on this measurement.
A: Approved. (Less common in daily tables than blank/NA, but sometimes used). It generally indicates the data has been reviewed and approved.
Other common symbols you might see:
E: Estimated (data was missing and filled in).
R: Revised (data was updated after initial publication).
Note: If the Value column is NA, it means no data was recorded for that day. If the Symbol column is NA, it means the recorded data is standard/good.
HYDAT daily data are in Cubic Meters per Second and Meters.
