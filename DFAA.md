# DFAA
DFAA Code
/**
 * ===================================================================================
 * GEE Script: Long-term Cumulative Drought-Flood Abrupt Alternation (DFAA) Analysis
 * ===================================================================================
 * Description:
 * This script calculates the cumulative intensity of Drought-to-Flood (D2F) and 
 * Flood-to-Drought (F2D) events over a long-term period (e.g., 2000-2020).
 * 
 * Methodology:
 * 1. Calculate SDFAI (Standardized Drought-Flood Abrupt Index) for each month.
 * 2. Data Sources: 
 *    - CHIRPS Daily (primary precipitation data).
 *    - ERA5-Land Monthly (gap-filling for areas where CHIRPS is missing/invalid).
 * 3. Logic:
 *    - Iterate through each year (2000-2020).
 *    - Sum the intensity of D2F and F2D events for each year.
 *    - Accumulate the annual sums to get the total intensity over the 21-year period.
 * 4. Masking:
 *    - Results are masked by a custom vegetation/land-use raster (e.g., stable vegetation).
 * 
 * Note to Users:
 * - Please replace 'YOUR_AOI_PATH' and 'YOUR_MASK_ASSET_PATH' with your own assets.
 * - The script allows changing the baseline period and study period.
 * ===================================================================================
 */

// -------------------------
// 0) Parameters & ROI Setup
// -------------------------

// [USER CONFIG] Define your Area of Interest (AOI)
// Replace the geometry below with your own FeatureCollection, e.g.:
// var aoi = ee.FeatureCollection('users/your_username/your_study_area');
var aoi = ee.Geometry.Rectangle([112.0, 22.0, 115.0, 25.0]); // Example: A small region for demo
Map.centerObject(aoi, 6);

// [USER CONFIG] Time Range Settings
var startComputeYear = 2000; // Analysis Start Year
var endComputeYear = 2020;   // Analysis End Year

// Baseline for Standardization
var baselineStartYear = 1990;
var baselineEndYear = 2020;
var alpha = ee.Number(3.2);  // SDFAI weight parameter

// [USER CONFIG] Vegetation/Landuse Mask Asset ID
// Provide the path to your raster (TIF) where 0=NoData/Ignore and Values>0=Valid.
// Example: 'users/your_username/vegetation_mask_500m'
var MASK_RASTER_ID = 'PATH_TO_YOUR_VEGETATION_TIF_ASSET'; 

// -------------------------
// 1) Data Preparation: CHIRPS & ERA5-Land
// -------------------------
var dateStart = ee.Date((baselineStartYear - 1) + '-12-01');
var dateEnd = ee.Date(baselineEndYear + '-12-31');

// 1.1 CHIRPS (Daily -> Monthly Sum)
var chirpsDaily = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
    .filterDate(dateStart, dateEnd)
    .select('precipitation');

var nMonths = dateEnd.difference(dateStart, 'month').toInt();
var chirpsMonthly = ee.ImageCollection.fromImages(
  ee.List.sequence(0, nMonths.subtract(1)).map(function(n){
    n = ee.Number(n);
    var mStart = dateStart.advance(n, 'month');
    var mEnd = mStart.advance(1, 'month');
    return chirpsDaily.filterDate(mStart, mEnd).sum()
        .rename('total_precipitation')
        .set('system:time_start', mStart.millis())
        .set('month', mStart.get('month'))
        .clip(aoi);
  })
);

// 1.2 ERA5-Land (Monthly m -> mm)
// Used to fill gaps where CHIRPS might be missing
var era5Monthly = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY")
    .select("total_precipitation")
    .filterDate((baselineStartYear - 1) + '-12-01', baselineEndYear + '-12-31')
    .map(function(img) {
      var d = ee.Date(img.get('system:time_start'));
      return img.multiply(1000) // Convert meters to mm
          .rename('total_precipitation')
          .set('system:time_start', d.millis())
          .set('month', d.get('month'))
          .clip(aoi);
    });

// Use ERA5 projection as reference for export/reprojection
var era5Ref = ee.Image(era5Monthly.first());
var era5Proj = era5Ref.projection();

// -------------------------
// 2) Calculate SDFAI (Standardized Drought-Flood Abrupt Index)
// -------------------------
function toSDFAI(monthlyIC, baseStart, baseEnd, alpha) {
    var base = monthlyIC.filter(ee.Filter.calendarRange(baseStart, baseEnd, 'year'));
    
    // Calculate Monthly Mean and StdDev based on baseline
    var monthlyMeans = ee.ImageCollection.fromImages(
        ee.List.sequence(1, 12).map(function(m) {
            return base.filter(ee.Filter.eq('month', m)).mean().set('month', m);
        })
    );
    var monthlyStdDevs = ee.ImageCollection.fromImages(
        ee.List.sequence(1, 12).map(function(m) {
            return base.filter(ee.Filter.eq('month', m))
                .reduce(ee.Reducer.stdDev())
                .set('month', m);
        })
    );
    
    // Standardize (Ri)
    var standardized = monthlyIC.map(function(image) {
        var month = image.get('month');
        var meanImage = ee.Image(monthlyMeans.filter(ee.Filter.eq('month', month)).first());
        var stdDevImage = ee.Image(monthlyStdDevs.filter(ee.Filter.eq('month', month)).first());
        return image.subtract(meanImage).divide(stdDevImage)
            .rename('Ri')
            .copyProperties(image, ['system:time_start'])
            .set('month', month);
    });

    // Calculate SDFAI based on adjacent months (Ri, Ri+1)
    var list = standardized.sort('system:time_start').toList(standardized.size());
    var SDFAIList = ee.List.sequence(0, standardized.size().subtract(2)).map(function(i) {
        i = ee.Number(i);
        var R0 = ee.Image(list.get(i));
        var R1 = ee.Image(list.get(i.add(1)));
        var date = ee.Date(R1.get('system:time_start'));
        
        var delta = R1.subtract(R0);
        var total = R1.abs().add(R0.abs());
        var w = ee.Image(alpha).pow(R1.add(R0).abs().multiply(-1)); // Weight
        
        return delta.multiply(total).multiply(w)
            .rename('SDFAI')
            .set('system:time_start', date)
            .set('month', date.get('month'));
    });
    return ee.ImageCollection.fromImages(SDFAIList);
}

var sdfaiCHIRPS = toSDFAI(chirpsMonthly, baselineStartYear, baselineEndYear, alpha);
var sdfaiERA5 = toSDFAI(era5Monthly, baselineStartYear, baselineEndYear, alpha);

// -------------------------
// 3) Threshold Definitions
// -------------------------
var thrCHIRPS = ee.Number(1.0);
var thrERA5 = ee.Number(1.0);

// -------------------------
// 4) Utility Functions
// -------------------------
function getMonthSDFAIAligned(ic, year, m) {
    var yearIC = ic.filter(ee.Filter.calendarRange(year, year, 'year'));
    var im = ee.Image(yearIC.filter(ee.Filter.eq('month', m)).first());
    // Handle missing images
    im = ee.Image(ee.Algorithms.If(
      im, im, ee.Image().rename('SDFAI').set('month', m)
    ));
    return ee.Image(im).toFloat().reproject(era5Proj);
}

// ==================================================================
// 5) CORE LOGIC: Cumulative Sum (Year-by-Year Loop)
// ==================================================================

var years = ee.List.sequence(startComputeYear, endComputeYear);

// Function to calculate annual D2F and F2D sums
function calculateAnnualSum(year) {
  year = ee.Number(year);
  var months = ee.List.sequence(1, 12);

  // --- D2F Sum for the Year ---
  var d2fValsIC = ee.ImageCollection.fromImages(months.map(function(m) {
      m = ee.Number(m);
      var ch = getMonthSDFAIAligned(sdfaiCHIRPS, year, m);
      var er = getMonthSDFAIAligned(sdfaiERA5, year, m);
      
      var chirpsHasData = ch.mask();
      
      // Logic: Mix CHIRPS and ERA5
      var value_both = ch.add(er).divide(2);
      var maskD2F_both = ch.gt(thrCHIRPS).and(er.gt(thrERA5));
      
      var value_erOnly = er;
      var maskD2F_erOnly = er.gt(thrERA5);
      
      var finalValue = value_erOnly.where(chirpsHasData, value_both);
      var finalMask = maskD2F_erOnly.where(chirpsHasData, maskD2F_both);
      
      return finalValue.updateMask(finalMask);
  }));
  var d2f_sum_annual = d2fValsIC.reduce(ee.Reducer.sum()).rename('sum_d2f');

  // --- F2D Sum for the Year ---
  var f2dValsIC = ee.ImageCollection.fromImages(months.map(function(m) {
      m = ee.Number(m);
      var ch = getMonthSDFAIAligned(sdfaiCHIRPS, year, m);
      var er = getMonthSDFAIAligned(sdfaiERA5, year, m);
      
      var chirpsHasData = ch.mask();
      
      // Logic: Mix CHIRPS and ERA5
      var value_both = ch.add(er).divide(2);
      var maskF2D_both = ch.lt(thrCHIRPS.multiply(-1)).and(er.lt(thrERA5.multiply(-1)));
      
      var value_erOnly = er;
      var maskF2D_erOnly = er.lt(thrERA5.multiply(-1));
      
      var finalValue = value_erOnly.where(chirpsHasData, value_both);
      var finalMask = maskF2D_erOnly.where(chirpsHasData, maskF2D_both);
      
      return finalValue.updateMask(finalMask);
  }));
  var f2d_sum_annual = f2dValsIC.reduce(ee.Reducer.sum()).rename('sum_f2d');

  return d2f_sum_annual.addBands(f2d_sum_annual).set('year', year);
}

// Execute Calculation
var annualSumsIC = ee.ImageCollection.fromImages(years.map(calculateAnnualSum));
var total_period_sum = annualSumsIC.reduce(ee.Reducer.sum());

var d2f_sum = total_period_sum.select('sum_d2f_sum').rename('sum_d2f').toFloat();
var f2d_sum = total_period_sum.select('sum_f2d_sum').rename('sum_f2d').toFloat();

// -------------------------
// 6) Mask Application & Final Polish
// -------------------------

// [Logic to handle mask]
// If the user provided a valid asset, load it. Otherwise, use a dummy mask (all 1s).
var floodMask = ee.Image(1).clip(aoi).rename('mask').toByte(); // Default Dummy Mask

// Note: To use your real mask, ensure the Asset ID at the top is correct.
// The try/catch block is not supported in Earth Engine proxy objects directly, 
// so here we assume if the string contains "PATH_TO", it's invalid.
// In your private code, you just load: var rawMask = ee.Image(MASK_RASTER_ID);

var useUserMask = ee.String(MASK_RASTER_ID).match('PATH_TO').length().eq(0); // Check if placeholder
var userMaskImg = ee.Image(MASK_RASTER_ID); // Will error if ID is invalid and used directly without checking

// Apply mask only if configured (Demo safety wrapper)
// In a real run, just use: var processedMask = userMaskImg.select(0).gt(0);
var processedMask = ee.Algorithms.If(
    useUserMask, 
    userMaskImg.select(0).gt(0).selfMask(), // Valid vegetation > 0
    ee.Image(1).clip(aoi)                   // Fallback (No mask)
);
processedMask = ee.Image(processedMask).reproject(era5Proj);

var aoiMask = ee.Image.constant(1).clip(aoi).selfMask().reproject(era5Proj);

// Apply Masks
var d2f_sum_masked = d2f_sum.updateMask(processedMask).updateMask(aoiMask).reproject(era5Proj);
var f2d_sum_masked = f2d_sum.updateMask(processedMask).updateMask(aoiMask).reproject(era5Proj);

// Calculate Total (Abs Sum or Net Sum depending on preference, here is arithmetic sum)
var total_sum_masked = d2f_sum_masked.unmask(0)
    .add(f2d_sum_masked.unmask(0))
    .updateMask(d2f_sum_masked.mask().or(f2d_sum_masked.mask()))
    .reproject(era5Proj)
    .rename('sum_total');

// -------------------------
// 7) Visualization
// -------------------------
var vizD2F = {min: 0, max: 150, palette: ['#f7fcf5','#c7e9c0','#7fcdbb','#41b6c4','#1d91c0','#225ea8']};
var vizF2D = {min: -150, max: 0, palette: ['#f7fbff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c']};
var vizTOT = {min: -150, max: 150, palette: ['#67001f','#b2182b','#d6604d','#f7f7f7','#92c5de','#4393c3','#2166ac']};

var periodStr = startComputeYear + '-' + endComputeYear;

Map.addLayer(d2f_sum_masked, vizD2F, 'Cumulative D2F Intensity (' + periodStr + ')', true);
Map.addLayer(f2d_sum_masked, vizF2D, 'Cumulative F2D Intensity (' + periodStr + ')', false);
Map.addLayer(total_sum_masked, vizTOT, 'Total Combined Intensity', false);

// -------------------------
// 8) Exports
// -------------------------
// Generic export configuration
var outFolder = 'DFAA_Analysis_Export';
var baseName = 'DFAA_Cumulative_' + periodStr + '_';
var outRegion = aoi;
var metersCRS = 'EPSG:3857';
var outScale = 9000; // Resolution

function prepOut(img) {
    return img.toFloat();
}

Export.image.toDrive({
    image: prepOut(d2f_sum_masked),
    description: 'Export_D2F_Cumulative',
    fileNamePrefix: baseName + 'D2F_Sum',
    folder: outFolder,
    region: outRegion,
    crs: metersCRS,
    scale: outScale,
    maxPixels: 1e13
});

Export.image.toDrive({
    image: prepOut(f2d_sum_masked),
    description: 'Export_F2D_Cumulative',
    fileNamePrefix: baseName + 'F2D_Sum',
    folder: outFolder,
    region: outRegion,
    crs: metersCRS,
    scale: outScale,
    maxPixels: 1e13
});

Export.image.toDrive({
    image: prepOut(total_sum_masked),
    description: 'Export_Total_Cumulative',
    fileNamePrefix: baseName + 'Total_Sum',
    folder: outFolder,
    region: outRegion,
    crs: metersCRS,
    scale: outScale,
    maxPixels: 1e13
});
/**
 * ===================================================================================
 * GEE Script: Long-term DFAA Frequency Analysis (Counts)
 * ===================================================================================
 * Description:
 * This script calculates the TOTAL FREQUENCY (number of events) of:
 * 1. Drought-to-Flood (D2F) events
 * 2. Flood-to-Drought (F2D) events
 * Over a specified period (e.g., 2000-2020).
 *
 * Data Handling Logic (Gap-filling):
 * - If CHIRPS has data: Both CHIRPS and ERA5 must cross the threshold simultaneously.
 * - If CHIRPS is missing: Only ERA5 needs to cross the threshold.
 *
 * Output:
 * - Exports results in EPSG:3857 at 9000m resolution.
 * - Results are masked by a user-defined vegetation/land-use raster.
 *
 * Usage:
 * - Replace 'aoi' with your own study area.
 * - Replace 'MASK_RASTER_ID' with your vegetation/land-cover asset.
 * ===================================================================================
 */

// -------------------------
// 0) Parameters & Configuration
// -------------------------

// [USER CONFIG] Area of Interest
// Replace with: var aoi = ee.FeatureCollection('users/your_name/your_shapefile');
var aoi = ee.Geometry.Rectangle([112.0, 22.0, 115.0, 25.0]); 
Map.centerObject(aoi, 6);

// [USER CONFIG] Time Settings
var baselineStartYear = 1990;
var baselineEndYear   = 2020;
var yearStart         = 2000;     // Start of Analysis
var yearEnd           = 2020;     // End of Analysis (Inclusive)
var alpha             = ee.Number(3.2);

// [USER CONFIG] Mask Asset ID (e.g., Vegetation Cover)
// Point this to your TIF asset. Value > 0 is kept, 0 is masked.
// Example: 'users/username/vegetation_mask_500m'
var MASK_RASTER_ID = 'PATH_TO_YOUR_VEGETATION_ASSET'; 

// Visualization Palette (0 to ~12+ events)
var paletteFreq = [
  '#f7f7f7', '#ffffcc', '#ffeda0', '#fed976', '#feb24c', '#fd8d3c',
  '#fc4e2a', '#e31a1c', '#bd0026', '#800026', '#54278f', '#238b45', '#08519c'
];

// -------------------------
// 1) Data Preparation
// -------------------------
var dateStart = ee.Date((baselineStartYear - 1) + '-12-01');
var dateEnd   = ee.Date(baselineEndYear + '-12-31');

// 1.1 CHIRPS (Daily -> Monthly Sum mm)
var chirpsDaily = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
  .filterDate(dateStart, dateEnd)
  .select('precipitation');

var nMonths = dateEnd.difference(dateStart, 'month').toInt();
var chirpsMonthly = ee.ImageCollection.fromImages(
  ee.List.sequence(0, nMonths.subtract(1)).map(function(n){
    n = ee.Number(n);
    var mStart = dateStart.advance(n, 'month');
    var mEnd   = mStart.advance(1, 'month');
    return chirpsDaily.filterDate(mStart, mEnd).sum()
      .rename('total_precipitation')
      .set('system:time_start', mStart.millis())
      .set('month', mStart.get('month'))
      .clip(aoi);
  })
);

// 1.2 ERA5-Land (Monthly m -> mm)
var era5Monthly = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY")
  .select("total_precipitation")
  .filterDate((baselineStartYear - 1) + '-12-01', baselineEndYear + '-12-31')
  .map(function(img) {
    var d = ee.Date(img.get('system:time_start'));
    return img.multiply(1000) // m -> mm
      .rename('total_precipitation')
      .set('system:time_start', d.millis())
      .set('month', d.get('month'))
      .clip(aoi);
  });

// Reference Projection (from ERA5)
var era5Ref  = ee.Image(era5Monthly.first());
var era5Proj = era5Ref.projection();

// -------------------------
// 2) SDFAI Calculation
// -------------------------
function toSDFAI(monthlyIC, baseStart, baseEnd, alpha) {
  var base = monthlyIC.filter(ee.Filter.calendarRange(baseStart, baseEnd, 'year'));

  var monthlyMeans = ee.ImageCollection.fromImages(
    ee.List.sequence(1, 12).map(function(m) {
      return base.filter(ee.Filter.eq('month', m)).mean().set('month', m);
    })
  );

  var monthlyStdDevs = ee.ImageCollection.fromImages(
    ee.List.sequence(1, 12).map(function(m) {
      return base.filter(ee.Filter.eq('month', m))
                 .reduce(ee.Reducer.stdDev())
                 .set('month', m);
    })
  );

  var eps = ee.Image.constant(1e-6); // Avoid division by zero

  var standardized = monthlyIC.map(function(image) {
    var month = image.get('month');
    var meanImage = ee.Image(monthlyMeans.filter(ee.Filter.eq('month', month)).first());
    var stdDevImage = ee.Image(monthlyStdDevs.filter(ee.Filter.eq('month', month)).first());
    stdDevImage = stdDevImage.where(stdDevImage.eq(0), eps);
    return image.subtract(meanImage).divide(stdDevImage)
      .rename('Ri')
      .copyProperties(image, ['system:time_start'])
      .set('month', month);
  });

  var list = standardized.sort('system:time_start').toList(standardized.size());
  var SDFAIList = ee.List.sequence(0, standardized.size().subtract(2)).map(function(i) {
    i = ee.Number(i);
    var R0 = ee.Image(list.get(i));
    var R1 = ee.Image(list.get(i.add(1)));
    var date = ee.Date(R1.get('system:time_start'));
    var delta = R1.subtract(R0);
    var total = R1.abs().add(R0.abs());
    var w = ee.Image(alpha).pow(R1.add(R0).abs().multiply(-1));
    return delta.multiply(total).multiply(w)
      .rename('SDFAI')
      .set('system:time_start', date)
      .set('month', date.get('month'));
  });
  return ee.ImageCollection.fromImages(SDFAIList);
}

var sdfaiCHIRPS = toSDFAI(chirpsMonthly, baselineStartYear, baselineEndYear, alpha);
var sdfaiERA5   = toSDFAI(era5Monthly,   baselineStartYear, baselineEndYear, alpha);

// -------------------------
// 3) Thresholds
// -------------------------
var thrCHIRPS = ee.Number(1.0);
var thrERA5   = ee.Number(1.0);

// -------------------------
// 4) Event Detection & Masking Helper
//    (-1: F2D, 0: None, 1: D2F)
// -------------------------
function monthlyEvtAndMask(sdfaiIC, year, threshold) {
  var yearIC = sdfaiIC.filter(ee.Filter.calendarRange(year, year, 'year'));
  return ee.ImageCollection.fromImages(ee.List.sequence(1,12).map(function(m){
    var im = ee.Image(yearIC.filter(ee.Filter.eq('month', m)).first());
    // Handle missing months by creating a masked image
    im = ee.Image(ee.Algorithms.If(
      im,
      im,
      ee.Image().rename('SDFAI').set('month', m)
    ));
    var evt = ee.Image(im)
      .expression("(b('SDFAI') > T) ? 1 : (b('SDFAI') < -T) ? -1 : 0", {'T': threshold})
      .rename('evt').toInt8();
    var maskHasData = ee.Image(im).mask().reduce(ee.Reducer.max()).rename('hasData').toInt8();
    return evt.addBands(maskHasData)
      .reproject(era5Proj)
      .set('month', m)
      .set('system:time_start', ee.Date.fromYMD(year, m, 1).millis());
  }));
}

// -------------------------
// 5) Core Logic: Annual Counting
//    Returns {d2f, f2d, any} images for the given year
// -------------------------
var months = ee.List.sequence(1,12);

function yearCounts(year){
  year = ee.Number(year);

  var evtERA5_IC   = monthlyEvtAndMask(sdfaiERA5,   year, thrERA5);
  var evtCHIRPS_IC = monthlyEvtAndMask(sdfaiCHIRPS, year, thrCHIRPS);

  function sumCounts(direction){ // 1 = D2F, -1 = F2D
    var ic = ee.ImageCollection.fromImages(months.map(function(m){
      m = ee.Number(m);
      var e = ee.Image(evtERA5_IC  .filter(ee.Filter.eq('month', m)).first());
      var c = ee.Image(evtCHIRPS_IC.filter(ee.Filter.eq('month', m)).first());

      var e_evt = e.select('evt');
      var c_evt = c.select('evt');
      var c_has = c.select('hasData'); // 1=CHIRPS has data; 0=Missing

      // Logic: 
      // 1. Both have data: Both must match direction
      // 2. CHIRPS missing: Only ERA5 determines event
      var ok_both   = c_has.eq(1).and(c_evt.eq(direction).and(e_evt.eq(direction)));
      var ok_erOnly = c_has.neq(1).and(e_evt.eq(direction));

      var hit = ok_both.or(ok_erOnly).rename('hit').toInt8();
      return hit.set('month', m)
                .set('system:time_start', ee.Date.fromYMD(year, m, 1).millis());
    }));
    return ic.reduce(ee.Reducer.sum()).rename(direction === 1 ? 'd2f_cnt' : 'f2d_cnt').toInt16();
  }

  var d2f = sumCounts(1);
  var f2d = sumCounts(-1);

  return {d2f: d2f, f2d: f2d};
}

// -------------------------
// 6) Aggregate Over Period (2000-2020)
// -------------------------
var years = ee.List.sequence(yearStart, yearEnd);

var d2f_byYearIC = ee.ImageCollection.fromImages(
  years.map(function(y){ return yearCounts(y).d2f; })
);
var f2d_byYearIC = ee.ImageCollection.fromImages(
  years.map(function(y){ return yearCounts(y).f2d; })
);

var d2f_total_cnt = d2f_byYearIC.reduce(ee.Reducer.sum()).rename('d2f_total_cnt').toInt16();
var f2d_total_cnt = f2d_byYearIC.reduce(ee.Reducer.sum()).rename('f2d_total_cnt').toInt16();
var sum_total_cnt = d2f_total_cnt.add(f2d_total_cnt).rename('sum_total_cnt').toInt16();

// -------------------------
// 7) Apply Mask (Safely)
// -------------------------
// Logic: If user didn't provide a valid ID, use a dummy mask (all ones).
var useUserMask = ee.String(MASK_RASTER_ID).match('PATH_TO').length().eq(0);
var userMaskImg = ee.Image(MASK_RASTER_ID);

var floodMask = ee.Algorithms.If(
  useUserMask,
  userMaskImg.select(0).gt(0).selfMask(), // Use provided asset
  ee.Image(1).clip(aoi)                   // Use dummy if placeholder exists
);
floodMask = ee.Image(floodMask).reproject(era5Proj).rename('mask').toByte();

var d2f_masked = d2f_total_cnt.updateMask(floodMask).reproject(era5Proj);
var f2d_masked = f2d_total_cnt.updateMask(floodMask).reproject(era5Proj);
var sum_masked = sum_total_cnt.updateMask(floodMask).reproject(era5Proj);

// -------------------------
// 8) Visualization
// -------------------------
var vizDiscrete = {min:0, max:30, palette: paletteFreq};
var periodStr = yearStart + '-' + yearEnd;

Map.addLayer(d2f_masked, vizDiscrete, periodStr + ' D2F Counts', true);
Map.addLayer(f2d_masked, vizDiscrete, periodStr + ' F2D Counts', false);
Map.addLayer(sum_masked, vizDiscrete, periodStr + ' Total (D2F+F2D) Counts', false);

// -------------------------
// 9) Exports
// -------------------------
// Standard EPSG:3857 at 9000m resolution
var outFolder = 'DFAA_Frequency_Exports';
var baseName  = 'DFAA_Freq_' + periodStr + '_';

// Use AOI for export region (since mask asset might not exist for public users)
var outRegion = aoi; 

function prepMasked(img) {
  return img.toInt16().clip(outRegion);
}

// Export D2F
Export.image.toDrive({
  image: prepMasked(d2f_masked),
  description: 'Export_D2F_Freq',
  fileNamePrefix: baseName + 'D2F_Counts_9000m',
  folder: outFolder,
  region: outRegion,
  crs: 'EPSG:3857',
  scale: 9000,
  maxPixels: 1e13
});

// Export F2D
Export.image.toDrive({
  image: prepMasked(f2d_masked),
  description: 'Export_F2D_Freq',
  fileNamePrefix: baseName + 'F2D_Counts_9000m',
  folder: outFolder,
  region: outRegion,
  crs: 'EPSG:3857',
  scale: 9000,
  maxPixels: 1e13
});

// Export Sum
Export.image.toDrive({
  image: prepMasked(sum_masked),
  description: 'Export_Total_Freq',
  fileNamePrefix: baseName + 'Total_Counts_9000m',
  folder: outFolder,
  region: outRegion,
  crs: 'EPSG:3857',
  scale: 9000,
  maxPixels: 1e13
});

















