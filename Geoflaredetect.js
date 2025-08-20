// ===============================================================================
// GeoFlare Detect: Aplikasi Analisis Kebakaran Hutan Interaktif di Rokan Hilir, Riau
// ===============================================================================

// --- 1. Koleksi Data dan Definisi Aset ---
var s2Sr = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var s2Clouds = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY');
var area = ee.FeatureCollection('projects/ee-fashihafirta1404/assets/Rokan_Hilir'); // Batas wilayah studi 
var firms = ee.ImageCollection('FIRMS'); // Data hotspot FIRMS (global, akan difilter nanti)

// Definisikan dataset validasi dan sampel pelatihan
var burned_2019 = ee.FeatureCollection('projects/ee-fashihafirta1404/assets/burned_2019');
var unburned_2019 = ee.FeatureCollection('projects/ee-fashihafirta1404/assets/unburned_2019');
var testvalidation19 = ee.FeatureCollection('projects/ee-fashihafirta1404/assets/Burnscar_2019_Clip');

var burned_2021 = ee.FeatureCollection('projects/ee-fashihafirta1404/assets/burned_2021');
var unburned_2021 = ee.FeatureCollection('projects/ee-fashihafirta1404/assets/unburned_2021');
var testvalidation21 = ee.FeatureCollection('projects/ee-fashihafirta1404/assets/Burnscar_2021_Clip');

var burned_2023 = ee.FeatureCollection('projects/ee-fashihafirta1404/assets/burned_2023');
var unburned_2023 = ee.FeatureCollection('projects/ee-fashihafirta1404/assets/unburned_2023');
var testvalidation23 = ee.FeatureCollection('projects/ee-fashihafirta1404/assets/Burnscar_2023_Clip');

var klasifikasi2019 = ee.Image('projects/ee-fashihafirta1404/assets/Klasifikasi_Kebakaran_2019');
var klasifikasi2021 = ee.Image('projects/ee-fashihafirta1404/assets/Klasifikasi_Kebakaran_2021');
var klasifikasi2023 = ee.Image('projects/ee-fashihafirta1404/assets/Klasifikasi_Kebakaran_2023');

// --- 2. Variabel Global dan Visualisasi ---u
var maxCloudProb = 25; // Probabilitas awan maksimum untuk masking

// SLD untuk dNBR
var sld_intervals =
  '<RasterSymbolizer>' +
    '<ColorMap type="intervals" extended="false">' +
      '<ColorMapEntry color="#ffffff" quantity="-500" label="Enhanced Regrowth High" />' +
      '<ColorMapEntry color="#7a8737" quantity="-250" label="Enhanced Regrowth Low" />' +
      '<ColorMapEntry color="#acbe4d" quantity="-100" label="Unburnt" />' +
      '<ColorMapEntry color="#0ae042" quantity="100" label="Low Severity" />' +
      '<ColorMapEntry color="#fff70b" quantity="270" label="Moderate-Low" />' +
      '<ColorMapEntry color="#ffaf38" quantity="440" label="Moderate-High" />' +
      '<ColorMapEntry color="#ff641b" quantity="660" label="High Severity" />' +
      '<ColorMapEntry color="#a41fd6" quantity="2000" label="No Data"/>' +
    '</ColorMap>' +
  '</RasterSymbolizer>';
  
var sld_dndvi =
  '<RasterSymbolizer>' +
    '<ColorMap type="intervals" extended="false">' +
      '<ColorMapEntry color="#008000" quantity="0.10" label="Tidak Berubah" />' +
      '<ColorMapEntry color="#99cc00" quantity="0.27" label="Keparahan Rendah" />' +
      '<ColorMapEntry color="#ffff00" quantity="0.44" label="Keparahan Sedang-Rendah" />' +
      '<ColorMapEntry color="#ff7f00" quantity="0.66" label="Keparahan Sedang-Tinggi" />' +
      '<ColorMapEntry color="#ff0000" quantity="1" label="Keparahan Tinggi" />' +
    '</ColorMap>' +
  '</RasterSymbolizer>';

// Palet warna untuk NDVI dan dNDVI
var ndviPalette = ['8c510a', 'd8b365', 'f6e8c3', 'c7eae5', '5ab4ac', '01665e'];
var dndviPalette = ['d9d9d9','ff0000', 'ff9900', 'ffff00', 'ffff66','ccff66', '99cc00', '339900', '006600'];

// Visualisasi klasifikasi area terbakar
var burnedAreaVis = {min: 0, max: 1, palette: ['#00ff00', '#FF0000']}; // Ungu: Tidak terbakar, Merah: Terbakar
var firmsVis = {
  min: 0,
  max: 100,
  palette: ['#ffff00', '#ff6600', '#ff0000']
};
var validationVis = {color: 'purple'}; // Warna untuk data validasi lapangan
var burnedVis = {color: 'red'};
var unburnedVis = {color: 'green'};

// --- 3. Fungsi-fungsi Utama ---
function maskClouds(img, maxCloudProb) {
  var cloudProb = ee.Image(img.get('cloud_mask')).select('probability');
  return img.updateMask(cloudProb.lt(maxCloudProb));
}

function maskEdges(img) {
  return img.updateMask(img.select('B8A').mask().updateMask(img.select('B9').mask()));
}

function maskWaterSCL(image) {
  // Pilih band SCL dari citra
  var scl = image.select('SCL');
  var maskThis = scl.eq(3) // Cloud Shadows
                   .or(scl.eq(6)); // Water
  var finalMask = maskThis.not(); 
  
  return image.updateMask(finalMask);
}

function processPrePostFire(tahun, area, maxCloudProb) {
  var preDate = ee.Date.fromYMD(ee.Number(tahun), 1, 1);
  var midDate = ee.Date.fromYMD(ee.Number(tahun), 6, 30);
  var postDate = ee.Date.fromYMD(ee.Number(tahun), 7, 1);
  var endDate = ee.Date.fromYMD(ee.Number(tahun), 12, 31);

  // Koleksi citra pre-fire
  var s2Pre = s2Sr.filterBounds(area).filterDate(preDate, midDate).map(maskEdges);
  var s2CloudsPre = s2Clouds.filterBounds(area).filterDate(preDate, midDate);
  var s2PreJoined = ee.Join.saveFirst('cloud_mask').apply({
    primary: s2Pre, secondary: s2CloudsPre, condition: ee.Filter.equals({ leftField: 'system:index', rightField: 'system:index' })
  });
  var s2PreMasked = ee.ImageCollection(s2PreJoined).map(function(img) {
    return maskClouds(img, maxCloudProb);
  }).map(maskWaterSCL).median().clip(area);

  // Koleksi citra post-fire
  var s2Post = s2Sr.filterBounds(area).filterDate(postDate, endDate).map(maskEdges);
  var s2CloudsPost = s2Clouds.filterBounds(area).filterDate(postDate, endDate);
  var s2PostJoined = ee.Join.saveFirst('cloud_mask').apply({
    primary: s2Post, secondary: s2CloudsPost, condition: ee.Filter.equals({ leftField: 'system:index', rightField: 'system:index' })
  });
  var s2PostMasked = ee.ImageCollection(s2PostJoined).map(function(img) {
    return maskClouds(img, maxCloudProb);
  }).map(maskWaterSCL).median().clip(area);

  // Perhitungan Indeks Vegetasi dan Kebakaran
  var ndviPre = s2PreMasked.normalizedDifference(['B8', 'B4']).rename('NDVI_pre');
  var ndviPost = s2PostMasked.normalizedDifference(['B8', 'B4']).rename('NDVI_post');
  var dndvi = ndviPre.subtract(ndviPost).rename('dNDVI');

  var nbrPre = s2PreMasked.normalizedDifference(['B8', 'B12']).multiply(1000).rename('NBR_pre');
  var nbrPost = s2PostMasked.normalizedDifference(['B8', 'B12']).multiply(1000).rename('NBR_post');
  var dnbr = nbrPre.subtract(nbrPost).rename('dNBR');

  var rdnbr = dnbr.divide(nbrPre.abs().add(1).sqrt()).rename('RdNBR');
  var rbr = dnbr.divide(nbrPre.add(1)).rename('RBR');

  var ndbi = s2PostMasked.normalizedDifference(['B11', 'B8']).rename('NDBI'); // Untuk potensi identifikasi permukiman

  return {
    preFire: s2PreMasked,
    postFire: s2PostMasked,
    ndviPre: ndviPre,
    ndviPost: ndviPost,
    dndvi: dndvi,
    nbrPre: nbrPre,
    nbrPost: nbrPost,
    dnbr: dnbr,
    rdnbr: rdnbr,
    rbr: rbr,
    ndbi: ndbi
  };
}

function trainAndClassify(year, postImage, dnbrImage, rdnbrImage, rbrImage, dndviImage, burnedSamples, unburnedSamples) {
  var featureImage = dnbrImage.addBands(rdnbrImage)
                               .addBands(rbrImage)
                               .addBands(dndviImage)
                               .addBands(postImage.select(['B4', 'B8', 'B11', 'B12']));

  var bands = ['dNBR', 'RdNBR', 'RBR', 'dNDVI',
              // 'NDVI_pre', 'NDVI_post',
              // 'NBR_pre', 'NBR_post', 'NDBI', 
              //'B2', 'B3', 
              'B4', 'B8', 'B11', 'B12'];

  // Mengambil sampel acak dari feature collection
  var sampledBurned = ee.FeatureCollection.randomPoints({
    region: burnedSamples.geometry(),
    points: 1000, // Jumlah titik sampel
    seed: 0
  }).map(function(f) { return f.set('class', 1); });

  var sampledUnburned = ee.FeatureCollection.randomPoints({
    region: unburnedSamples.geometry(),
    points: 1000,
    seed: 0
  }).map(function(f) { return f.set('class', 0); });

  var sampel = sampledBurned.merge(sampledUnburned);

  // Mengambil sampel citra untuk pelatihan
  var training = featureImage.sampleRegions({
    collection: sampel, properties: ['class'], scale: 20
  });

  // Membagi data menjadi data pelatihan dan pengujian
  var withRandom = training.randomColumn('random');
  var split = 0.8;
  var trainData = withRandom.filter(ee.Filter.lt('random', split));
  var testData = withRandom.filter(ee.Filter.gte('random', split));

  // Melatih classifier Random Forest
  var classifier = ee.Classifier.smileRandomForest({
    numberOfTrees: 100,
    variablesPerSplit: null,
    minLeafPopulation: 5,
    bagFraction: 0.5,
    maxNodes: null,
    seed: 0
  }).train({
    features: trainData,
    classProperty: 'class',
    inputProperties: bands
  });

  // Mengklasifikasikan citra
  var classified = featureImage.classify(classifier);

  // Penilaian Akurasi Pelatihan
  var trainAccuracy = classifier.confusionMatrix();

  // Penilaian Akurasi Pengujian
  var testAccuracy = testData.classify(classifier).errorMatrix('class', 'classification');

  return {
    classified: classified,
    classifier: classifier,
    trainAccuracy: trainAccuracy,
    testAccuracy: testAccuracy,
    featureImage: featureImage,
    bands: bands
  };
}

// FUNGSI LEGENDA SKALA WARNA
function makeColorRampLegend(colors, labels, title) {
  var legendPanel = ui.Panel({
    style: {
      padding: '0px',
      margin: '0px',
      shown: false // Sembunyikan secara default
    }
  });
  legendPanel.add(ui.Label(title, {fontWeight: 'bold'}));
  var legendContent = ui.Panel({ layout: ui.Panel.Layout.Flow('vertical') });

  for (var i = 0; i < colors.length; i++) {
    var colorBox = ui.Label({
      style: { backgroundColor: colors[i], padding: '8px', margin: '0px 4px 0px 0px' }
    });
    var description = ui.Label({
      value: labels[i], style: { margin: '0px 0px 4px 0px' }
    });
    legendContent.add(ui.Panel([colorBox, description], ui.Panel.Layout.Flow('horizontal')));
  }
  legendPanel.add(legendContent);
  return legendPanel;
}

// FUNGSI LEGENDA KATEGORI
function makeCategoricalLegend(colors, labels, title) {
  var legendPanel = ui.Panel({
    style: {
      padding: '0px',
      margin: '0px',
      shown: false // Sembunyikan secara default
    }
  });
  legendPanel.add(ui.Label(title, { fontWeight: 'bold' }));
  var legendContent = ui.Panel({ layout: ui.Panel.Layout.Flow('vertical') });

  for (var i = 0; i < colors.length; i++) {
    var colorBox = ui.Label({
      style: { backgroundColor: colors[i], padding: '8px', margin: '0px 4px 0px 0px' }
    });
    var description = ui.Label({
      value: labels[i], style: { margin: '0px 0px 4px 0px' }
    });
    legendContent.add(ui.Panel([colorBox, description], ui.Panel.Layout.Flow('horizontal')));
  }
  legendPanel.add(legendContent);
  return legendPanel;
}

// 1. Legenda Batas Admin
// Membuat legenda batas admin
var adminLegend = null; 

// Membuat legenda batas admin
function makeAdminLegend() {
  // Jika legenda sudah ada, jangan buat lagi.
  if (adminLegend) {
    return adminLegend; 
  }

  var panel = ui.Panel({
    style: {
      position: 'bottom-right',
      padding: '8px 15px',
      stretch: 'horizontal'
      //backgroundColor: 'rgba(255, 255, 255, 0.8)'
    }
  });

  panel.add(ui.Label('Batas Administrasi', { fontWeight: 'bold' }));

  var legendContent = ui.Panel({ 
    layout: ui.Panel.Layout.Flow('horizontal'), 
    style: {stretch: 'horizontal'} 
  });

  var colorBox = ui.Label({
    style: {
      backgroundColor: 'black',
      padding: '8px',
      //margin: '0px 4px 0px 0px'
    }
  });

  var description = ui.Label('Wilayah Kajian');

  legendContent.add(ui.Panel([colorBox, description], ui.Panel.Layout.Flow('horizontal')));
  panel.add(legendContent);

  // Simpan panel legenda yang baru dibuat ke variabel global
  //adminLegend = panel; 
  
  return panel;
}

// 2. Legenda Titik Panas FIRMS 
function makeFirmsLegend(year) {
  // Gunakan color ramp legend yang lebih informatif
  var legendPanel = ui.Panel({
    style: {
      position: 'bottom-right', // Letak di peta
      padding: '8px 15px',
      shown:true,
      backgroundColor: 'rgba(255, 255, 255, 0.8)'
    }
  });

  // Tambahkan judul legenda
  var legendTitle = ui.Label('Tingkat Kepercayaan Titik Panas ' + year, {
    fontWeight: 'bold',
    margin: '0 0 4px 0'
  });
  legendPanel.add(legendTitle);

  // Warna dan label sesuai skala confidence FIRMS
  var colors = ['#ffff00', '#ff6600', '#ff0000'];
  var labels = ['Rendah (0% - 30%)', 'Sedang (30% - 80%)', 'Tinggi (80% - 100%)'];

  // Tambahkan isi legenda
  for (var i = 0; i < colors.length; i++) {
    var colorBox = ui.Label('', {
      backgroundColor: colors[i],
      padding: '8px',
      margin: '0px 4px 0px 0px'
    });
    var desc = ui.Label(labels[i], { margin: '0px 0px 4px 0px' });
    legendPanel.add(ui.Panel([colorBox, desc], ui.Panel.Layout.Flow('horizontal')));
  }
  return legendPanel;
}

// 3. Legenda Data Validasi Lapangan
function makeValidationLegend(year) {
  var legendPanel = ui.Panel({
    style: {
      position: 'bottom-right', // Posisi di dalam peta
      padding: '8px 15px',
      shown: true,
      backgroundColor: 'rgba(255, 255, 255, 0.8)'
    }
  });
  legendPanel.add(ui.Label('Data Validasi Lapangan ' + year, { fontWeight: 'bold' }));
  var legendContent = ui.Panel({ layout: ui.Panel.Layout.Flow('horizontal') });
  var colorBox = ui.Label({
    style: { backgroundColor: validationVis.color, padding: '8px', margin: '0px 4px 0px 0px' }
  });
  var description = ui.Label({
    value: 'Area Terbakar Terverifikasi', style: { margin: '0px 0px 4px 0px' }
  });
  legendContent.add(ui.Panel([colorBox, description], ui.Panel.Layout.Flow('horizontal')));
  legendPanel.add(legendContent);
  return legendPanel;
}

// 4. Legenda Burned Sampel (Area Terbakar)
function makeBurnedLegend(year) {
  // Gunakan color ramp legend yang lebih informatif
  var legendPanel = ui.Panel({
    style: {
      position: 'bottom-right', // Letak di peta
      padding: '8px 15px',
      shown:true,
      //backgroundColor: 'rgba(255, 255, 255, 0.8)'
    }
  });
  legendPanel.add(ui.Label(' Sampel Area Terbakar ' + year, { fontWeight: 'bold' }));
  var legendContent = ui.Panel({ layout: ui.Panel.Layout.Flow('horizontal') });
  var colorBox = ui.Label({
    style: { backgroundColor: burnedVis.color, padding: '8px', margin: '0px 4px 0px 0px' }
  });
  var description = ui.Label({
    value: 'Sampel Area Terbakar', style: { margin: '0px 0px 4px 0px' }
  });
  legendContent.add(ui.Panel([colorBox, description], ui.Panel.Layout.Flow('horizontal')));
  legendPanel.add(legendContent);
  return legendPanel;
}

// 4. Legenda Unburned Sampel (Area Tidak Terbakar)
function makeUnburnedLegend(year) {
  // Gunakan color ramp legend yang lebih informatif
  var legendPanel = ui.Panel({
    style: {
      position: 'bottom-right', // Letak di peta
      padding: '8px 15px',
      shown:true,
      //backgroundColor: 'rgba(255, 255, 255, 0.8)'
    }
  });
  legendPanel.add(ui.Label(' Sampel Area Tidak Terbakar ' + year, { fontWeight: 'bold' }));
  var legendContent = ui.Panel({ layout: ui.Panel.Layout.Flow('horizontal') });
  var colorBox = ui.Label({
    style: { backgroundColor: unburnedVis.color, padding: '8px', margin: '0px 4px 0px 0px' }
  });
  var description = ui.Label({
    value: 'Sampel Area Tidak Terbakar', style: { margin: '0px 0px 4px 0px' }
  });
  legendContent.add(ui.Panel([colorBox, description], ui.Panel.Layout.Flow('horizontal')));
  legendPanel.add(legendContent);
  return legendPanel;
}

// 5. Legenda Klasifikasi Kebakaran
function makeKlasifikasiLegend(year) {
  // Gunakan color ramp legend yang lebih informatif
  var legendPanel = ui.Panel({
    style: {
      position: 'bottom-right', // Letak di peta
      padding: '8px 15px',
      shown:true
    }
  });
  // Tambahkan judul legenda
  var legendTitle = ui.Label(' Klasifikasi Kebakaran ' + year, {
    fontWeight: 'bold',
    margin: '0 0 4px 0'
  });
  legendPanel.add(legendTitle);
  var colors = ['#00ff00', '#FF0000'];
  var labels = ['Area Tidak Terbakar', 'Area Terbakar'];
  // Tambahkan isi legenda
  for (var i = 0; i < colors.length; i++) {
    var colorBox = ui.Label('', {
      backgroundColor: colors[i],
      padding: '8px',
      margin: '0px 4px 0px 0px'
    });
    var desc = ui.Label(labels[i], { margin: '0px 0px 4px 0px' });
    legendPanel.add(ui.Panel([colorBox, desc], ui.Panel.Layout.Flow('horizontal')));
  }
  return legendPanel;
}
 
// 6. Legenda dNBR & RdNBR & RBR
function makedNBRLegend(year) {
  // Gunakan color ramp legend yang lebih informatif
  var legendPanel = ui.Panel({
    style: {
      position: 'bottom-right', // Letak di peta
      padding: '8px 15px',
      shown:true
    }
  });

  // Tambahkan judul legenda
  var legendTitle = ui.Label('Kelas dNBR ' + year, {
    fontWeight: 'bold',
    margin: '0 0 4px 0'
  });
  legendPanel.add(legendTitle);

  // Daftar warna dan label sesuai sld_intervals
  var entries = [
    {color: '#7a8737', label: 'Pertumbuhan Kembali Tinggi'},
    {color: '#acbe4d', label: 'Pertumbuhan Kembali Rendah'},
    {color: '#0ae042', label: 'Tidak Terbakar'},
    {color: '#fff70b', label: 'Kebakaran Ringan'},
    {color: '#ffaf38', label: 'Kebakaran Sedang-Rendah'},
    {color: '#ff641b', label: 'Kebakaran Sedang-Tinggi'},
    {color: '#a41fd6', label: 'Kebakaran Tinggi'},
    {color: '#ffffff', label: 'Tidak Ada Data'}
  ];

  entries.forEach(function(entry) {
    var colorBox = ui.Label('', {
      backgroundColor: entry.color,
      padding: '8px',
      margin: '0 4px 4px 0'
    });
    var label = ui.Label(entry.label, {margin: '0 0 4px 0'});
    legendPanel.add(ui.Panel([colorBox, label], ui.Panel.Layout.Flow('horizontal')));
  });

  return legendPanel;
}

// 7. Legenda RdNBR 
function makerdnbrLegend(year) {
  // Gunakan color ramp legend yang lebih informatif
  var legendPanel = ui.Panel({
    style: {
      position: 'bottom-right', // Letak di peta
      padding: '8px 15px',
      shown:true
    }
  });

  // Tambahkan judul legenda
  var legendTitle = ui.Label('Kelas RdNBR ' + year, {
    fontWeight: 'bold',
    margin: '0 0 4px 0'
  });
  legendPanel.add(legendTitle);

  // Daftar warna dan label sesuai sld_intervals
  var entries = [
    {color: '#7a8737', label: 'Pertumbuhan Kembali Tinggi'},
    {color: '#acbe4d', label: 'Pertumbuhan Kembali Rendah'},
    {color: '#0ae042', label: 'Tidak Terbakar'},
    {color: '#fff70b', label: 'Kebakaran Ringan'},
    {color: '#ffaf38', label: 'Kebakaran Sedang-Rendah'},
    {color: '#ff641b', label: 'Kebakaran Sedang-Tinggi'},
    {color: '#a41fd6', label: 'Kebakaran Tinggi'},
    {color: '#ffffff', label: 'Tidak Ada Data'}
  ];

  entries.forEach(function(entry) {
    var colorBox = ui.Label('', {
      backgroundColor: entry.color,
      padding: '8px',
      margin: '0 4px 4px 0'
    });
    var label = ui.Label(entry.label, {margin: '0 0 4px 0'});
    legendPanel.add(ui.Panel([colorBox, label], ui.Panel.Layout.Flow('horizontal')));
  });

  return legendPanel;
}

// 8. Legenda RBR
function makerbrLegend(year) {
  // Gunakan color ramp legend yang lebih informatif
  var legendPanel = ui.Panel({
    style: {
      position: 'bottom-right', // Letak di peta
      padding: '8px 15px',
      shown:true
    }
  });

  // Tambahkan judul legenda
  var legendTitle = ui.Label('Kelas RBR ' + year, {
    fontWeight: 'bold',
    margin: '0 0 4px 0'
  });
  legendPanel.add(legendTitle);

  // Daftar warna dan label sesuai sld_intervals
  var entries = [
    {color: '#7a8737', label: 'Pertumbuhan Kembali Tinggi'},
    {color: '#acbe4d', label: 'Pertumbuhan Kembali Rendah'},
    {color: '#0ae042', label: 'Tidak Terbakar'},
    {color: '#fff70b', label: 'Kebakaran Ringan'},
    {color: '#ffaf38', label: 'Kebakaran Sedang-Rendah'},
    {color: '#ff641b', label: 'Kebakaran Sedang-Tinggi'},
    {color: '#a41fd6', label: 'Kebakaran Tinggi'},
    {color: '#ffffff', label: 'Tidak Ada Data'}
  ];

  entries.forEach(function(entry) {
    var colorBox = ui.Label('', {
      backgroundColor: entry.color,
      padding: '8px',
      margin: '0 4px 4px 0'
    });
    var label = ui.Label(entry.label, {margin: '0 0 4px 0'});
    legendPanel.add(ui.Panel([colorBox, label], ui.Panel.Layout.Flow('horizontal')));
  });

  return legendPanel;
}

// 9. Legenda dNDVI
function makedNDVILegend(year) {
  // Gunakan color ramp legend yang lebih informatif
  var legendPanel = ui.Panel({
    style: {
      position: 'bottom-right', // Letak di peta
      padding: '8px 15px',
      shown:true
    }
  });

  // Tambahkan judul legenda
  var legendTitle = ui.Label('Kelas dNDVI ' + year, {
    fontWeight: 'bold',
    margin: '0 0 4px 0'
  });
  legendPanel.add(legendTitle);

  // Daftar warna dan label sesuai sld_intervals
  var entries = [
  {color: '#008000', label: 'Tidak Berubah'},
  {color: '#99cc00', label: 'Keparahan Rendah'},
  {color: '#ffff00', label: 'Keparahan Sedang-Rendah'},
  {color: '#ff7f00', label: 'Keparahan Sedang-Tinggi'},
  {color: '#ff0000', label: 'Keparahan Tinggi'}
  ];

  entries.forEach(function(entry) {
    var colorBox = ui.Label('', {
      backgroundColor: entry.color,
      padding: '8px',
      margin: '0 4px 4px 0'
    });
    var label = ui.Label(entry.label, {margin: '0 0 4px 0'});
    legendPanel.add(ui.Panel([colorBox, label], ui.Panel.Layout.Flow('horizontal')));
  });

  return legendPanel;
}

// --- 4. Proses Awal untuk Data Grafik (sekali saat aplikasi dimuat) ---

var allYearsData = {};
var years = [2019, 2021, 2023];

years.forEach(function(year) {
  var burnedSamples, unburnedSamples;
  // Perbaiki logika untuk memilih sampel pelatihan yang benar
  if (year === 2019) { burnedSamples = burned_2019; unburnedSamples = unburned_2019; }
  else if (year === 2021) { burnedSamples = burned_2021; unburnedSamples = unburned_2021; }
  else if (year === 2023) { burnedSamples = burned_2023; unburnedSamples = unburned_2023; }

  var processed = processPrePostFire(year, area, maxCloudProb);
  var classification = trainAndClassify(
    year, processed.postFire, processed.dnbr, processed.rdnbr, processed.rbr, processed.dndvi,
    // processed.ndviPre, processed.ndviPost, processed.nbrPre, processed.nbrPost, processed.ndbi,
    burnedSamples, unburnedSamples
  );
  allYearsData[year.toString()] = {
    processed: processed,
    classification: classification,
    validation: (year === 2019) ? testvalidation19 : (year === 2021) ? testvalidation21 : testvalidation23
  };
});

// --- 5. Implementasi UI (Disesuaikan dengan Kerangka Anda) ---
ui.root.clear(); // Bersihkan semua elemen root

// Inisialisasi Peta
var uiMap = ui.Map();
uiMap.setOptions("ROADMAP");
uiMap.setCenter(100.789, 2.121, 9); // Fokus ke Rokan Hilir

// Panel legenda (ditambahkan ke uiMap langsung)
var legendPanel = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px',
    shown: false,
    backgroundColor: 'rgba(255, 255, 255, 0.8)'
  }
});
uiMap.add(legendPanel); // Tambahkan legenda ke peta

// Panel Kontrol Kiri
var controlPanel = ui.Panel({
  layout: ui.Panel.Layout.Flow('vertical'),
  style: {
    width: '380px',
    padding: '15px',
    border: '1px solid #ccc',
    stretch: 'horizontal'
  }
});

// Panel Hasil Kanan
var resultsPanel = ui.Panel({
  layout: ui.Panel.Layout.Flow('vertical'),
  style: {
    width: '320px',
    padding: '15px',
    border: '1px solid #ccc',
    shown: false,
    stretch: 'horizontal'
  }
});

// Tambahkan panel ke root
ui.root.insert(0, controlPanel);
ui.root.insert(1, uiMap);
ui.root.insert(2, resultsPanel);


/// PANEL KIRI KONTROL
// Judul dan Deskripsi
controlPanel.add(ui.Label('GeoFlare Detect', {   fontWeight: 'bold',
  fontSize: '24px',
  color: '#8C0505', //
  textAlign: 'center',
  stretch: 'horizontal',
  padding: '8px 0',
  margin: '10px 0 5px 0'
  
}));
controlPanel.add(ui.Label(
  'Aplikasi yang dirancang untuk visualisasi dan analisis spasial kebakaran hutan dan lahan secara interaktif di wilayah Kabupaten Rokan Hilir, Provinsi Riau. Aplikasi ini berasal dari implementasi machine learning dengan algoritma Random Forest yang digunakan untuk mendeteksi dan memetakan area terbakar selama periode tahun 2019, 2021 dan 2023, data hotspot dari FIRMS, hasil klasifiikasi area terbakar berdasarkan citra Sentinel-2, serta indeks diferensial dari vegetasi dan kebakaran seperti NDVI dan NBR. Aplikasi ini memungkinkan pengguna dalam memantau kondisi kebakaran secara historis dan spasial serta mendukung proses pengambilan keputusan berbasis data.',
  { fontSize: '12px', margin: '10px 0px', textAlign: 'justify' }
));
//controlPanel.add(ui.Label('---', {}, null, {margin: '10px 0'}));


// BUTTON: Instruksi Penggunaan
var instructionsButton = ui.Button('Instruksi Penggunaan Earth Engine Apps', null, false, {
  margin: '5px 0px',
  stretch: 'horizontal'
});

var instructionsPanel = ui.Panel({style: {shown: false, fontSize: '12px'}});

// Judul Video Tutorial
instructionsPanel.add(ui.Label({
  value: '🎥 Video Tutorial:',
  style: {fontWeight: 'bold', fontSize: '12px', margin: '6px 0 4px 0'}
}));

// Link Video
instructionsPanel.add(ui.Label({
  value: 'Tonton di YouTube',
  style: {color: 'blue',fontSize: '12px', margin: '2px 0 2px 6px'},
  targetUrl: 'https://youtu.be/jcMBp1whnao' // Ganti dengan link YouTube Anda
}));

// Judul Tata Cara Penggunaan
instructionsPanel.add(ui.Label({
  value: '❓ Tata Cara Penggunaan:',
  style: {fontWeight: 'bold', fontSize: '12px', margin: '6px 0 4px 0'}
}));

// Langkah-langkah Penggunaan
var steps = [
  '1. Pilih tahun analisis dari dropdown.',
  '2. Ceklis layer peta yang ingin ditampilkan.',
  '3. Hasil analisis tampil di panel kanan (Evaluasi Hasil).',
  '4. Harap tunggu proses pemuatan peta, komputasi berjalan sesuai kecepatan jaringan internet.',
  '5. Gunakan zoom/pan untuk eksplorasi peta.',
  '6. Link ke hasil klasifikasi peta tersedia di informasi Glosarium.'
];

steps.forEach(function(step) {
  instructionsPanel.add(ui.Label({
    value: step,
    style: {fontSize: '12px', margin: '2px 0 2px 6px'}
  }));
});

instructionsButton.onClick(function() {
  instructionsPanel.style().set('shown', !instructionsPanel.style().get('shown'));
});
controlPanel.add(instructionsButton);
controlPanel.add(instructionsPanel);

// Analisis Kebakaran per Tahun
// controlPanel.add(ui.Label('Analisis Kebakaran per Tahun:', {fontWeight: 'bold'}));
var yearSelector = ui.Select({
  items: ['2019', '2021', '2023'],
  placeholder: 'Analisis Kebakaran per Tahun',
  onChange: function(year) {
    displayYearAnalysis(parseInt(year));
  },
  style: {stretch: 'horizontal', margin: '5px 0px'}
});
controlPanel.add(yearSelector);

// Layer Peta
//controlPanel.add(ui.Label('Layer Peta:', {fontWeight: 'bold'}));
var layerToggleButton = ui.Button({
  label: 'Layer Peta Deteksi Kebakaran Hutan dan Lahan',
  style: {stretch: 'horizontal', margin: '5px 0px'},
  onClick: function() {
    var isShown = layerOptionsPanel.style().get('shown');
    layerOptionsPanel.style().set('shown', !isShown);
    layerToggleButton.setLabel(isShown ? 'Layer Peta Deteksi Kebakaran Hutan dan Lahan' : 'Sembunyikan Layer');
  }
});

var layerOptionsPanel = ui.Panel({
  layout: ui.Panel.Layout.Flow('vertical'),
  style: {shown: false}  // default sembunyi
});

controlPanel.add(layerToggleButton);
controlPanel.add(layerOptionsPanel);

var layerOptionsPanel = ui.Panel({
  layout: ui.Panel.Layout.Flow('vertical'),
  style: {shown: false}
});
controlPanel.add(layerOptionsPanel); // Tambahkan ke panel utama

layerOptionsPanel.add(ui.Label('Pilih Layer Peta:', {fontWeight: 'bold', margin: '10px 0px', textAlign: 'justify'}));

// Checkbox Batas Administrasi
layerOptionsPanel.add(ui.Label('Batas Administrasi:', {fontWeight: 'bold', margin: '5px 0'}));
var adminBoundaryCheckbox = ui.Checkbox({
  label: 'Batas Administrasi'
});
layerOptionsPanel.add(adminBoundaryCheckbox);
//adminLayerPanel.add(adminBoundaryCheckbox); // atau controlPanel.add(...)

// Group Citra
layerOptionsPanel.add(ui.Label('Citra:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
var preFireRgbCheckbox = ui.Checkbox('Pre-Fire RGB', false);
var postFireRgbCheckbox = ui.Checkbox('Post-Fire RGB', false);
layerOptionsPanel.add(preFireRgbCheckbox);
layerOptionsPanel.add(postFireRgbCheckbox);

// Group Indeks Spektral
layerOptionsPanel.add(ui.Label('Indeks Spektral:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
var dnbrCheckbox = ui.Checkbox('dNBR', false);
var rdnbrCheckbox = ui.Checkbox('RdNBR', false);
var rbrCheckbox = ui.Checkbox('RBR', false);
var dndviCheckbox = ui.Checkbox('dNDVI', false);

layerOptionsPanel.add(dnbrCheckbox);
layerOptionsPanel.add(rdnbrCheckbox);
layerOptionsPanel.add(rbrCheckbox);
layerOptionsPanel.add(dndviCheckbox);

// Group untuk dNBR, RdNBR, RBR
var dNBRLayerPanel = ui.Panel();
layerOptionsPanel.add(dNBRLayerPanel);

var rdnbrLayerPanel = ui.Panel();
layerOptionsPanel.add(rdnbrLayerPanel);

var rbrLayerPanel = ui.Panel();
layerOptionsPanel.add(rbrLayerPanel);

var dNDVILayerPanel = ui.Panel();
layerOptionsPanel.add(dNDVILayerPanel);

// Group untuk Titik Panas Firms (akan diisi dinamis)
layerOptionsPanel.add(ui.Label('Titik Panas (Hotspot) FIRMS:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
var firmsLayerPanel = ui.Panel();
layerOptionsPanel.add(firmsLayerPanel);

// Group untuk Data Validasi Lapangan (akan diisi dinamis)
layerOptionsPanel.add(ui.Label('Data Validasi Lapangan:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
var validationLayerPanel = ui.Panel();
layerOptionsPanel.add(validationLayerPanel);

// Group untuk Sampel Lapangan (akan diisi dinamis)
layerOptionsPanel.add(ui.Label(' Sampel:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
var burnedLayerPanel = ui.Panel();
layerOptionsPanel.add(burnedLayerPanel);
var unburnedLayerPanel = ui.Panel();
layerOptionsPanel.add(unburnedLayerPanel);

// Group untuk Klasifikasi Kebakaran (akan diisi dinamis)
layerOptionsPanel.add(ui.Label('Hasil Klasifikasi Kebakaran:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
var klasifikasiLayerPanel = ui.Panel();
layerOptionsPanel.add(klasifikasiLayerPanel);

// Informasi Pengolahan Data
//controlPanel.add(ui.Label('Informasi Pengolahan Data: ', {fontWeight: 'bold'}));

// Panel isi glosarium (default tersembunyi)
var informasiPanel = ui.Panel({style: {shown: false, fontSize: '12px'}});

var informasiButton = ui.Button({
  label: 'Informasi Pengolahan Data',
  style: {margin: '5px 0px', stretch: 'horizontal'},
  onClick: function() {
    var isShown = informasiPanel.style().get('shown');
    informasiPanel.style().set('shown', !isShown);
    informasiButton.setLabel(isShown ? ' Informasi Pengolahan Data' : 'Sembunyikan Informasi Pengolahan Data');
  }
});

// Panel Informasi

// Bagian 1: Data yang Digunakan
informasiPanel.add(ui.Label({
  value: '📘 Data yang Digunakan:',
  style: {fontWeight: 'bold', fontSize: '12px', margin: '6px 0 4px 0'}
}));

var dataList = [
  '1. Citra Sentinel-2 (Level-2A) tahun 2019, 2021, dan 2023',
  '2. Data Titik Panas (Hotspot) Kebakaran FIRMS ',
  '3. Data batas administrasi Kabupaten Rokan Hilir (Sumber: BIG)',
  '4. Data Penutup Lahan Kabupaten Rokan Hilir (Sumber: BIG)',
  '5. Data Area Bekas Kebakaran Hutan dan Lahan di Indonesia (Sumber: KLHK)'
];

dataList.forEach(function(item) {
  informasiPanel.add(ui.Label({
    value: item,
    style: {fontSize: '12px', margin: '2px 0 2px 6px', textAlign: 'justify'}
  }));
});

// Bagian 2: Metode yang Digunakan
informasiPanel.add(ui.Label({
  value: '🛠️ Metode yang Digunakan:',
  style: {fontWeight: 'bold', fontSize: '12px', margin: '6px 0 4px 0'}
}));

var metodeList = [
  '1. Algoritma Utama: Menggunakan Random Forest untuk memodelkan klasifikasi area terbakar.',
  '2. Sumber Data: Memanfaatkan citra satelit Sentinel-2 (pre & post) dan indeks turunan seperti dNBR dan dNDVI untuk mendeteksi area terbakar.',
  '3. Persiapan Data: Melakukan proses masking awan pada citra, serta menggunakan data batas administrasi dan area bekas kebakaran sebagai pendukung.',
  '4. Pelatihan dan Uji Model: Model dilatih menggunakan sampel yang dibagi menjadi 80% data pelatihan dan 20% data pengujian.',
  '5. Evaluasi: Kinerja model dinilai menggunakan Confusion Matrix yang menghasilkan metrik seperti Overall Accuracy, Kappa Coefficient, Producer Accuracy, User Accuracy.'
];

metodeList.forEach(function(item) {
  informasiPanel.add(ui.Label({
    value: item,
    style: {fontSize: '12px', margin: '2px 0 2px 6px', textAlign: 'justify'}
  }));
});

// Bagian 3: Glosarium
informasiPanel.add(ui.Label({
  value: '📘 Glosarium:',
  style: {fontWeight: 'bold', fontSize: '12px', margin: '6px 0 4px 0'}
}));

var glosariumList = [
  '● Sentinel-2: Satelit penginderaan jauh yang menyediakan citra resolusi tinggi, sering digunakan untuk pemantauan lingkungan.',
  '● Pre dan Post Citra: Citra satelit yang diambil sebelum (pre) dan sesudah (post) sebuah kejadian, seperti kebakaran, untuk menganalisis perubahannya.',
  '● FIRMS: Sistem yang menyediakan data titik panas (hotspot) dan kebakaran aktif dari satelit seperti MODIS dan VIIRS secara real-time.',
  '● Random Forest: Sebuah model machine learning yang menggunakan banyak pohon keputusan untuk klasifikasi yang stabil dan akurat.',
  '● Sekumpulan kecil data yang diambil dari populasi atau keseluruhan data untuk mewakilinya dalam analisis',
  '● NDVI (Normalized Difference Vegetation Index): Indeks yang mengukur kesehatan dan kerapatan vegetasi, digunakan untuk membedakan vegetasi dari area terbakar.',
  '● NBR (Normalized Burn Ratio): Indeks yang dihitung dari citra satelit untuk mengukur tingkat keparahan kebakaran dan kondisi vegetasi.',
  '● dNBR (differenced NBR):  Nilai yang dihitung dari selisih NBR post-fire dan pre-fire untuk mengukur tingkat keparahan kebakaran secara lebih akurat.',
  '● RdNBR (Relative dNBR): dNBR yang dinormalisasi antar wilayah.',
  '● RBR (Relative Burn Ratio): Variasi lain untuk keparahan.',
  '● Training Data: Data yang digunakan untuk melatih model agar belajar mengenali pola.',
  '● Testing Data: Data terpisah yang digunakan untuk menguji performa akhir model setelah pelatihan.',
  '● Confusion Matrix: Tabel yang merangkum hasil prediksi model, menunjukkan prediksi yang benar dan salah untuk setiap kelas.',
  '● Overall Accuracy: Persentase total piksel yang diklasifikasikan dengan benar.',
  '● Kappa Coefficient: Metrik akurasi yang lebih kuat, memperhitungkan prediksi yang benar secara kebetulan.',
  '● Producer Accuracy: Mengukur seberapa baik model mendeteksi piksel dari suatu kelas yang sebenarnya (Recall).',
  '● User Accuracy: Mengukur seberapa sering prediksi model untuk suatu kelas benar-benar akurat (Precision).'
];

glosariumList.forEach(function(item) {
  informasiPanel.add(ui.Label({
    value: item,
    style: {fontSize: '12px', margin: '2px 0 2px 6px', textAlign: 'justify'}
  }));
});

// Tambahan: Link Download
informasiPanel.add(ui.Label({
  value: '[Klik: Download Peta Klasifikasi Kebakaran Model Random Forest]',
  style: {fontSize: '12px', color: 'blue', textAlign: 'left', margin: '2px 0 2px 6px'},
  targetUrl: 'https://drive.google.com/drive/folders/1mo-fI_yl8ILfsdT4Oku-kcPfVcoQbW9S?usp=sharing'
}));

// Tambahkan tombol dan panel ke controlPanel
controlPanel.add(informasiButton);
controlPanel.add(informasiPanel);

// Garis horizontal lurus (seperti <hr>)
var horizontalLine = ui.Panel({
  style: {
    border: '0.5px solid #999', // warna dan ketebalan garis
    margin: '10px 0px',         // jarak atas/bawah dan kiri/kanan
    height: '1px',               // tinggi garis
     stretch: 'horizontal'
  }
});

controlPanel.add(horizontalLine);

// Credit
controlPanel.add(ui.Label('Credit', {fontWeight: 'bold'}));
controlPanel.add(ui.Label('● Aplikasi ini dibuat oleh Fashiha Firta Prakasa sebagai syarat Proyek Akhir (PA), serta dibimbing oleh  Dr. Nur Mohammad Farda, S.Si., M.Cs.', {fontSize: '12px', textAlign: 'justify' }));
controlPanel.add(ui.Label('● Apabila terdapat pertanyaan atau saran mengenai platform ini silakan menghubungi:', {fontSize: '12px', textAlign: 'justify' }));
controlPanel.add(ui.Label({
  value: 'fashihaafirtaa@gmail.com',
  style: {fontSize: '12px', color: 'blue', margin: '2px 0 2px 6px', textAlign: 'center'},
  targetUrl: 'mailto:fashihaafirtaa@gmail.com'
}));
controlPanel.add(ui.Label({
  value: 'LinkedIn: linkedin.com/in/fashihafirta',
  style: {fontSize: '12px', color: 'blue', margin: '2px 0 2px 6px', textAlign: 'center'},
  targetUrl: 'https://www.linkedin.com/in/fashihafirta'  // ganti dengan profil kamu
}));
controlPanel.add(ui.Label(
  'Program Studi Sarjana Terapan Sistem Informasi Geografis\n' +
  'Departemen Teknologi Kebumian\n' +
  'Sekolah Vokasi\n' +
  'Universitas Gadjah Mada\n' +
  '2025',
  {
    fontSize: '12px',
    fontWeight: 'bold',
    textAlign: 'center',
    stretch: 'horizontal',
    whiteSpace: 'pre-line',
    margin: '10px 0 0 0'
  }
));

/// PANEL KANAN HASIL
// Konten dinamis akan ditambahkan di sini oleh displayYearAnalysis

/// FUNGSI UNTUK MENGATUR LAYER DAN HASIL
var currentMapLegends = []; // Untuk melacak legenda yang ditambahkan langsung ke Map

function clearMapLayers() {
  var layersToRemove = [];
  uiMap.layers().forEach(function(layer) {
    if (layer.getName() !== 'Batas Administrasi') { // Jangan hapus batas administrasi
      layersToRemove.push(layer);
    }
  });
  layersToRemove.forEach(function(layer) {
    uiMap.remove(layer);
  });
  currentMapLegends = [];
}

function clearYearResults() {
  resultsPanel.clear(); // Bersihkan hasil di panel kanan
  //resultsPanel.add(ui.Label('Hasil Analisis Kebakaran', { fontWeight: 'bold', fontSize: '18px', textAlign: 'center', stretch: 'horizontal' }));
  //resultsPanel.add(ui.Label('---', {}, null, {margin: '10px 0'}));
}
 // Memperbarui status checkbox dan layer peta.

function addLayerAndBindCheckbox(image, visParams, name, checkbox, legendCreator) {
  // Pastikan legendCreator tidak undefined jika tidak diberikan
  legendCreator = legendCreator || function() { return ui.Panel({style: {shown: false}}); };

  var layerName = name + ' ' + yearSelector.getValue();
  // Set initial layer state based on checkbox value
  if (checkbox.getValue()) {
      uiMap.addLayer(image, visParams, layerName);
      if (legendCreator) {
          var newLegend = legendCreator(yearSelector.getValue());
          //uiMap.add(newLegend);
          currentMapLegends.push(newLegend);
      }
  }

  // Bind onChange event for the checkbox
  checkbox.onChange(function(checked) {
    var currentYear = yearSelector.getValue();
    var specificLayerName = name + ' ' + currentYear;
    
        var existingLayer = null;
    // Iterasi melalui layer yang ada di peta untuk menemukan yang cocok
    for (var i = 0; i < uiMap.layers().length(); i++) {
        var layer = uiMap.layers().get(i);
        if (layer.getName() === specificLayerName) {
            existingLayer = layer;
            break; // Berhenti jika sudah ditemukan
        }
    }
    
    if (checked) {
      if (!existingLayer) { // Add only if not already present
        uiMap.addLayer(image, visParams, specificLayerName);
      } else {
        existingLayer.setShown(true);
      }
      if (legendCreator) {
      // Hapus semua legenda lama terlebih dahulu
      currentMapLegends.forEach(function(lg) {
        uiMap.remove(lg);
      });
      currentMapLegends = [];

      // Tambah legenda baru
      var newLegend = legendCreator(currentYear);
      uiMap.add(newLegend);
      currentMapLegends.push(newLegend);
    }
    } else {
      if (existingLayer) {
        uiMap.remove(existingLayer);
      }

      // Remove associated legend if it exists and is for this layer
      currentMapLegends = currentMapLegends.filter(function(legend) {
      var firstWidget = legend.widgets().get(0);
      if (firstWidget && typeof firstWidget.getValue === 'function') {
          var legendTitle = firstWidget.getValue();
          if (legendTitle && legendTitle.toString().indexOf(name) !== -1) {
          uiMap.remove(legend);
          return false;
          }
      }
          return true;
      });
    }
  });
}

function updateLayerVisibility(year) {
  var data = allYearsData[year.toString()];
  if (!data) return;

  // Hapus semua layer terkait tahun sebelumnya dan legenda
  clearMapLayers();
  
addLayerAndBindCheckbox(
  area.style({color: 'black', fillColor: '00000000'}),
  {},  // visParams kosong karena sudah styled
  'Batas Administrasi',
  adminBoundaryCheckbox,
  makeAdminLegend // bisa null jika tak mau ada legend
);


  // Reset all relevant checkboxes to false first (except admin boundary)
  [preFireRgbCheckbox, postFireRgbCheckbox, 
  dnbrCheckbox,rdnbrCheckbox, rbrCheckbox, 
   dndviCheckbox, 
   ]
   .forEach(function(cb) {
     cb.setValue(false, false); // Set value to false tanpa memicu onChange
   });


  // Re-create and bind FIRMS and Validation checkboxes for the current year
  firmsLayerPanel.clear();

var firmsCheckbox = ui.Checkbox({
  label: 'Titik Panas FIRMS (' + year + ')',
  value: false
});
firmsLayerPanel.add(firmsCheckbox);


// Layer & Legend Titik Panas FIRMS
addLayerAndBindCheckbox(
  firms.filterDate(year + '-01-01', year + '-12-31').filterBounds(area).select('confidence').median().clip(area),
  firmsVis,
  'Tingkat Kepercayaan Titik Panas',
  firmsCheckbox,
  makeFirmsLegend //
);

// Layer Data Validasi
  validationLayerPanel.clear();
  var validationDataset;
  if (year === 2019) validationDataset = testvalidation19;
  else if (year === 2021) validationDataset = testvalidation21;
  else if (year === 2023) validationDataset = testvalidation23;

  var validationCheckbox = ui.Checkbox({
    label: 'Data Validasi Lapangan (' + year + ')',
    value: false
  });
  validationLayerPanel.add(validationCheckbox);
  addLayerAndBindCheckbox(validationDataset, validationVis, 'Data Validasi Lapangan', validationCheckbox, makeValidationLegend);

// Layer Sampel Area Terbakar
burnedLayerPanel.clear();
  var burnedDataset;
  if (year === 2019) burnedDataset = burned_2019;
  else if (year === 2021) burnedDataset = burned_2021;
  else if (year === 2023) burnedDataset = burned_2023;

  var burnedCheckbox = ui.Checkbox({
    label: 'Sampel Area Terbakar (' + year + ')',
    value: false
  });
  burnedLayerPanel.add(burnedCheckbox);
  addLayerAndBindCheckbox(burnedDataset, burnedVis, 'Sampel Area Terbakar', burnedCheckbox, makeBurnedLegend);
  
  // Layer Sampel Area Tidak Terbakar
unburnedLayerPanel.clear();
  var unburnedDataset;
  if (year === 2019) unburnedDataset = unburned_2019;
  else if (year === 2021) unburnedDataset = unburned_2021;
  else if (year === 2023) unburnedDataset = unburned_2023;

  var unburnedCheckbox = ui.Checkbox({
    label: 'Sampel Area Tidak Terbakar (' + year + ')',
    value: false
  });
  unburnedLayerPanel.add(unburnedCheckbox);
  addLayerAndBindCheckbox(unburnedDataset, unburnedVis, 'Sampel Area Tidak Terbakar', unburnedCheckbox, makeUnburnedLegend);

// Layer Sampel Area Terbakar
klasifikasiLayerPanel.clear();
  var klasifikasiDataset;
  if (year === 2019) klasifikasiDataset = klasifikasi2019;
  else if (year === 2021) klasifikasiDataset = klasifikasi2021;
  else if (year === 2023) klasifikasiDataset = klasifikasi2023;

  var klasifikasiCheckbox = ui.Checkbox({
    label: 'Klasifikasi Kebakaran (' + year + ')',
    value: false
  });
  klasifikasiLayerPanel.add(klasifikasiCheckbox);
  addLayerAndBindCheckbox(klasifikasiDataset, burnedAreaVis, 'Klasifikasi Kebakaran', klasifikasiCheckbox, makeKlasifikasiLegend);

// Layer citra & indeks 
  preFireRgbCheckbox.setLabel('Pre-Fire RGB (' + year + ')');
  postFireRgbCheckbox.setLabel('Post-Fire RGB (' + year + ')');
  dNBRLayerPanel.clear();
  dnbrCheckbox.setLabel('dNBR (' + year + ')');
  rdnbrLayerPanel.clear();
  rdnbrCheckbox.setLabel('RdNBR (' + year + ')');
  rbrLayerPanel.clear();
  rbrCheckbox.setLabel('RBR (' + year + ')');
  dNDVILayerPanel.clear();
  dndviCheckbox.setLabel('dNDVI (' + year + ')');
  //classifiedBurnedAreaCheckbox.setLabel('Klasifikasi Kebakaran (' + year + ')');

  // Bind other layer checkboxes
  addLayerAndBindCheckbox(data.processed.preFire, {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000, gamma: 1.2}, 'Pre-Fire RGB', preFireRgbCheckbox);
  addLayerAndBindCheckbox(data.processed.postFire, {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000, gamma: 1.2}, 'Post-Fire RGB', postFireRgbCheckbox);
  addLayerAndBindCheckbox(data.processed.dnbr.sldStyle(sld_intervals), {}, 'dNBR', dnbrCheckbox, makedNBRLegend);
  addLayerAndBindCheckbox(data.processed.rdnbr.sldStyle(sld_intervals), {}, 'RdNBR', rdnbrCheckbox, makerdnbrLegend);
  addLayerAndBindCheckbox(data.processed.rbr.sldStyle(sld_intervals), {}, 'RBR', rbrCheckbox, makerbrLegend);
  addLayerAndBindCheckbox(data.processed.dndvi.sldStyle(sld_dndvi), {}, 'dNDVI', dndviCheckbox, makedNDVILegend);
  //addLayerAndBindCheckbox(data.classification.classified, burnedAreaVis, 'Klasifikasi Kebakaran', classifiedBurnedAreaCheckbox);
}

// --- FUNGSI UNTUK MENGATUR LAYER DAN HASIL ---

// ===============================================================================
// Fungsi untuk menghitung total hotspot
// ===============================================================================

function countTotalHotspot(firmsCollection, year, geometry) {
  // Filter koleksi FIRMS berdasarkan tahun
  var yearlyFirms = firmsCollection.filterDate(year + '-01-01', year + '-12-31');

  // Gabungkan semua citra menjadi satu citra mosaik
  var combined = yearlyFirms.map(function(img) {
    // Pastikan band 'confidence' ada dan ambil hanya piksel dengan confidence > 0
    return img.select('confidence').updateMask(img.select('confidence').gt(0));
  }).mosaic();

  // Hitung jumlah piksel terdeteksi sebagai hotspot
  var count = combined.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: geometry, 
    scale: 1000,
    maxPixels: 1e13
  }).get('confidence');

  return ee.Feature(null, {year: year, count: count});
}

// ===============================================================================
// Fungsi Tampilkan Grafik Luas Terbakar (Direvisi untuk menambahkan ke panel spesifik)
// ===============================================================================
function tampilkanGrafikTitikPanas(targetPanelForChart) {
  var resolusiKlasifikasi = 20;

var hotspot2019 = ee.Number(countTotalHotspot(firms, 2019, Riau).get('count'));
  var hotspot2021 = ee.Number(countTotalHotspot(firms, 2021, Riau).get('count'));
  var hotspot2023 = ee.Number(countTotalHotspot(firms, 2023, Riau).get('count'));

  // 2. Hitung luas terbakar validasi dalam hektar
  var luasValidasi2019 = testvalidation19.geometry().area().divide(10000); // ha
  var luasValidasi2021 = testvalidation21.geometry().area().divide(10000); // ha
  var luasValidasi2023 = testvalidation23.geometry().area().divide(10000); // ha

  // Buat FeatureCollection untuk ComboChart
  var fiturTren = ee.FeatureCollection([
    ee.Feature(null, {
      tahun: 2019,
      hotspot: hotspot2019,
      luas_validasi: luasValidasi2019
    }),
    ee.Feature(null, {
      tahun: 2021,
      hotspot: hotspot2021,
      luas_validasi: luasValidasi2021
    }),
    ee.Feature(null, {
      tahun: 2023,
      hotspot: hotspot2023,
      luas_validasi: luasValidasi2023
    })
  ]);

  // Buat ComboChart (Hotspot sebagai garis, Luas Terbakar Validasi sebagai batang)
  var chartHotspot = ui.Chart.feature.byFeature(fiturTren, 'tahun', ['luas_validasi', 'hotspot'])
    .setChartType('ComboChart')
    .setOptions({
      title: 'Perbandingan Titik Panas dan Luas Validasi Area Terbakar',
      hAxis: {title: 'Tahun'},
      vAxes: {
        0: {title: 'Luas Validasi (ha)'},
        1: {title: 'Jumlah Titik Panas'}
      },
      series: {
        0: {type: 'bars', targetAxisIndex: 0, color: '#f57c00'}, // Luas validasi
        1: {type: 'line', targetAxisIndex: 1, color: '#1976d2', lineWidth: 3, pointSize: 5} // Titik panas
      },
      legend: {position: 'top'},
      height: 320,
      width: 420
    });
  targetPanelForChart.add(chartHotspot);
}

// Fungsi untuk menampilkan grafik batang luas kebakaran
function tampilkanGrafikLuasTerbakar(targetPanelForChart) {
  // Data luas kebakaran per tahun (ha)
  var tahun = ['2019', '2021', '2023'];
  var luasKebakaran = [0.8, 0.32, 0.68]; // Ganti dengan data aktualmu

  // Buat FeatureCollection dari data
  var dataChart = ee.FeatureCollection(
    tahun.map(function(t, i){
      return ee.Feature(null, {
        'Tahun': t,
        'LuasKebakaran': luasKebakaran[i]
      });
    })
  );

  // Buat chart batang
  var chartLuas = ui.Chart.feature.byFeature({
    features: dataChart,
    xProperty: 'Tahun',
    yProperties: ['LuasKebakaran']
  })
  .setChartType('ColumnChart')
  .setOptions({
    title: 'Luas Kebakaran per Tahun (ha)',
    hAxis: {title: 'Tahun'},
    vAxis: {title: 'Luas (ha)'},
    legend: {position: 'none'},
    colors: ['#d73027']
  });

  // Tampilkan chart ke dalam panel
  targetPanelForChart.add(chartLuas);
}

// Fungsi utama untuk menampilkan chart di panel target
function displayComparisonChart(targetPanel, callback) {
  targetPanel.add(ui.Label('Grafik Perbandingan Luas Terbakar per Tahun:', {
    fontWeight: 'bold',
    margin: '10px 0 5px 0'
  }));
  
  tampilkanGrafikLuasTerbakar(targetPanel); // Tambahkan grafik ke panel

  callback(); // Callback bisa langsung dipanggil karena proses sinkron
}

// ===============================================================================
// Fungsi Utama: displayYearAnalysis
// ===============================================================================
function displayYearAnalysis(currentYear) {
  resultsPanel.style().set('shown', true);
  clearYearResults();
  updateLayerVisibility(currentYear); // Pastikan ini berfungsi dengan 'currentYear'

if (adminBoundaryCheckbox.getValue()) {
  // Pastikan 'adminBoundary' didefinisikan untuk tahun yang dipilih atau secara global
  uiMap.addLayer(adminBoundary, {color: 'black'}, 'Batas Administrasi');
}
if (preFireRgbCheckbox.getValue()) {
  // Pastikan 'preFireImage' didefinisikan dan memiliki visParams yang benar
  uiMap.addLayer(preFireImage, preFireRgbVis, 'Pre-Fire RGB ' + year);
}
  var data = allYearsData[currentYear.toString()];
  if (!data) {
    resultsPanel.add(ui.Label('Data untuk tahun ' + currentYear + ' tidak tersedia.'));
    return;
  }

  resultsPanel.add(ui.Label(' Evaluasi Tahun ' + currentYear + '', {fontWeight: 'bold', fontSize: '18px', textAlign: 'center', stretch: 'horizontal' }));

  // --- Helper Function untuk Membuat Bagian Hasil yang Dapat Di-toggle ---
  function createToggleableSection(title, contentFunction) {
    var contentPanel = ui.Panel([], ui.Panel.Layout.Flow('vertical'));
    contentPanel.style().set('shown', false); // Awalnya tersembunyi
    var contentLoaded = false; // Flag untuk melacak apakah konten sudah dimuat
    var buttonLabel = '' + title;

    var button = ui.Button({
      label: buttonLabel,
      onClick: function() {
        if (!contentPanel.style().get('shown')) { // Jika tersembunyi
          if (!contentLoaded) { // Jika konten belum dimuat
            contentPanel.clear(); // Bersihkan pesan "Memuat..." jika ada
            //contentPanel.add(ui.Label('Memuat ' + title + '...', {fontStyle: 'italic', fontSize: '11px', color: 'grey'}));
            
            contentFunction(contentPanel, function() { // Panggil fungsi pemuatan konten
              contentPanel.style().set('shown', true); // Tampilkan panel konten
              button.setLabel('Sembunyikan ' + title);
              contentLoaded = true;
            });
          } else { // Konten sudah dimuat, tinggal tampilkan
            contentPanel.style().set('shown', true);
            button.setLabel('Sembunyikan ' + title);
          }
        } else { // Jika sudah terlihat
          contentPanel.style().set('shown', false); // Sembunyikan panel konten
          button.setLabel('' + title);
        }
      },
      style: {stretch: 'horizontal', margin: '5px 0'}
    });

    resultsPanel.add(button);
    resultsPanel.add(contentPanel); // Tambahkan panel konten (kosong dan tersembunyi)
  }

  // =============================================================================
  // Definisi Fungsi Konten untuk Setiap Bagian (Menambahkan ke 'targetPanel' yang diberikan)
  // =============================================================================

  // --- Akurasi Training Model ---
  function displayTrainAccuracy(targetPanel, callback) {
    targetPanel.add(ui.Label(' Akurasi Training Model:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
    var dataTraining = {
      2019: 0.9678362573099415,
      2021: 0.9802569727358195,
      2023: 0.9714786089567176
    };
    var dataKappa = {
    2019: 0.9263885261652215,
    2021: 0.9469853737932911,
    2023: 0.9224787183872182
    };
      var confusionMatrices = {
    2019: ee.Array([[1598,29], [76,734]]),
    2021: ee.Array([[2333,33], [48,755]]),
    2023: ee.Array([[2955,45], [93,904]])
    };
    
    var prodUserAccuracy = {
    2019: ee.Array([
      [0.9821, 0.9545],     // Kelas 0 (Tidak Terbakar)
      [0.9061, 0.9619]     // Kelas 1 (Terbakar)
    ]),
    2021: ee.Array([
      [0.9860, 0.9798],
      [0.9402, 0.9581]
    ]),
    2023: ee.Array([
      [0.985, 0.9694],     
      [0.9067,  0.9525]
    ])
};

    // Ambil nilai berdasarkan currentYear
    var akurasi = dataTraining[currentYear];
    var kappa = dataKappa[currentYear];
    var matrix = confusionMatrices[currentYear];
    var prodUserDict = prodUserAccuracy[currentYear];
    
    var formattedAkurasi = (akurasi !== null && akurasi !== undefined) ? akurasi.toFixed(4) : 'N/A';
    var formattedAkurasiPercent = (akurasi !== null && akurasi !== undefined) ? (akurasi * 100).toFixed(2) + '%' : 'N/A';
    var formattedKappa = (kappa !== null && kappa !== undefined) ? kappa.toFixed(4) : 'N/A';
    //var formattedProdUserPercent = (prodUserDict !== null && prodUserDict !== undefined) ? (prodUserDict * 100).toFixed(2) + '%' : 'N/A';
    
    targetPanel.add(ui.Label(' • Overall Accuracy ' +  ': ' + formattedAkurasi+ ' (' + formattedAkurasiPercent + ')'));
    targetPanel.add(ui.Label(' • Kappa Coefficient   ' +  ': ' + formattedKappa ));
    targetPanel.add(ui.Label('Confusion Matrix:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
    
    function createConfMatrixChart(matrix, currentYear) {
    return ui.Chart.array.values({
        array: matrix,
        axis: 0,
        xLabels: ['Aktual: Tidak Terbakar (0)', 'Aktual: Terbakar (1)']
      })
      .setSeriesNames(['Prediksi: Tidak Terbakar (0)', 'Prediksi: Terbakar (1)'])
      .setChartType('Table')
      .setOptions({
        title: 'Confusion Matrix (Training ' + currentYear + ')',
        hAxis: {title: 'Prediksi'},
        vAxis: {title: 'Aktual'},
        colors: ['#cccccc', '#cccccc'],
        legend: {position: 'none'},
        fontSize: 10,
        height: 200,
        width: 350
      });
  }
  
  // Tambahkan confusion matrix hanya untuk tahun yang dipilih
  targetPanel.add(createConfMatrixChart(matrix, currentYear));

  targetPanel.add(ui.Label('Akurasi Spesifik:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
    
function createProdUserChart(prodUserDict, currentYear) {
  // Ambil array dan ubah ke persen
  var accuracyArray = prodUserDict[currentYear];
  var accuracyPercent = accuracyArray.multiply(100);

  // Buat label untuk menjelaskan isi tabel
  var explanation = ui.Label({
    value: 'Kolom 0: Producer Accuracy (%)\nKolom 1: User Accuracy (%)',
    style: { fontSize: '11px', color: 'gray', margin: '0 0 4px 0' }
  });

  // Buat chart tipe tabel
  var chart = ui.Chart.array.values({
      array: accuracyPercent,
      axis: 0,
      xLabels: ['Tidak Terbakar (0)', 'Terbakar (1)']
    })
    .setSeriesNames(['Akurasi Producer', 'Akurasi User'])
    .setChartType('Table')
    .setOptions({
      title: 'Producer & User Accuracy (Training) - ' + currentYear,
      fontSize: 12,
      legend: { position: 'none' },
      width: 360,
      height: 200
    });

  // Gabungkan label dan chart dalam satu panel
  return ui.Panel([explanation, chart]);
}

// Contoh pemanggilan (pastikan 'produs' dan 'currentYear' sudah didefinisikan):
targetPanel.add(createProdUserChart(prodUserAccuracy, currentYear));
  
  callback();
}

  // --- Akurasi Testing Model ---
  function displayTestAccuracy(targetPanel, callback) {
    targetPanel.add(ui.Label('Akurasi Testing Model:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
    var dataTesting = {
      2019: 0.9438943894389439,
      2021: 0.9666254635352287,
      2023: 0.9803680981595092
    };
    var dataKappa = {
    2019: 0.875,
    2021: 0.9104974820630284,
    2023: 0.9462666886434812
    };
     var confusionMatricesTesting = {
    2019: ee.Array([[365,8], [22,168]]),
    2021: ee.Array([[582,12], [16,181]]),
    2023: ee.Array([[596,7], [19,181]])
    };
    var prodUserAccuracyTesting = {
    2019: ee.Array([
      [0.9785, 0.9431],     // Kelas 0 (Tidak Terbakar)
      [0.8842, 0.9545]     // Kelas 1 (Terbakar)
    ]),
    2021: ee.Array([
      [0.9797, 0.9732],
      [0.9187, 0.9378]
    ]),
    2023: ee.Array([
      [0.9883, 0.9691],     
      [0.905, 0.9627]
    ])
};
    // Ambil nilai berdasarkan currentYear
    var akurasi = dataTesting[currentYear];
    var kappa = dataKappa[currentYear];
    var matrix = confusionMatricesTesting[currentYear];
    var prodUserDict = prodUserAccuracyTesting[currentYear];
    
    var formattedAkurasi = (akurasi !== null && akurasi !== undefined) ? akurasi.toFixed(4) : 'N/A';
    var formattedAkurasiPercent = (akurasi !== null && akurasi !== undefined) ? (akurasi * 100).toFixed(2) + '%' : 'N/A';
    var formattedKappa = (kappa !== null && kappa !== undefined) ? kappa.toFixed(4) : 'N/A';
    
    targetPanel.add(ui.Label(' • Overall Accuracy ' +  ': ' + formattedAkurasi+ ' (' + formattedAkurasiPercent + ')'));
    targetPanel.add(ui.Label(' • Kappa Coefficient  ' +  ': ' + formattedKappa ));
    targetPanel.add(ui.Label('Confusion Matrix:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
    
    function createConfMatrixChartTesting(matrix, currentYear) {
    return ui.Chart.array.values({
        array: matrix,
        axis: 0,
        xLabels: ['Aktual: Tidak Terbakar (0)', 'Aktual: Terbakar (1)']
      })
      .setSeriesNames(['Prediksi: Tidak Terbakar (0)', 'Prediksi: Terbakar (1)'])
      .setChartType('Table')
      .setOptions({
        title: 'Confusion Matrix (Training ' + currentYear + ')',
        hAxis: {title: 'Prediksi'},
        vAxis: {title: 'Aktual'},
        colors: ['#cccccc', '#cccccc'],
        legend: {position: 'none'},
        fontSize: 10,
        height: 200,
        width: 350
      });
  }
  
  // Tambahkan confusion matrix hanya untuk tahun yang dipilih
  targetPanel.add(createConfMatrixChartTesting(matrix, currentYear));

  targetPanel.add(ui.Label('Akurasi Spesifik:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
    
function createProdUserChartTest(prodUserDict, currentYear) {
  // Ambil array dan ubah ke persen
  var accuracyArray = prodUserDict[currentYear];
  var accuracyPercent = accuracyArray.multiply(100);

  // Buat label untuk menjelaskan isi tabel
  var explanation = ui.Label({
    value: 'Kolom 0: Producer Accuracy (%)\nKolom 1: User Accuracy (%)',
    style: { fontSize: '11px', color: 'gray', margin: '0 0 4px 0' }
  });

  // Buat chart tipe tabel
  var chart = ui.Chart.array.values({
      array: accuracyPercent,
      axis: 0,
      xLabels: ['Tidak Terbakar (0)', 'Terbakar (1)']
    })
    .setSeriesNames(['Akurasi Producer', 'Akurasi User'])
    .setChartType('Table')
    .setOptions({
      title: 'Producer & User Accuracy (Training) - ' + currentYear,
      fontSize: 12,
      legend: { position: 'none' },
      width: 360,
      height: 200
    });

  // Gabungkan label dan chart dalam satu panel
  return ui.Panel([explanation, chart]);
  }

// Contoh pemanggilan (pastikan 'produs' dan 'currentYear' sudah didefinisikan):
  targetPanel.add(createProdUserChartTest(prodUserAccuracyTesting, currentYear));
    callback();
  }

  // --- Feature Importance ---
function displayFeatureImportance(targetPanel, callback) {
  targetPanel.add(ui.Label('Feature Importance:', {
    fontWeight: 'bold',
    margin: '10px 0 5px 0'
  }));

  // Data importance final per tahun
  var importanceData = {
    2019: {
      B11: 8.860564789654513,
      B12: 11.07128884217745,
      B4: 6.385015547448218,
      B8: 10.292536559654241,
      RBR: 8.218259928465248,
      RdNBR: 7.115953290608159,
      dNBR: 4.9558573730647595,
      dNDVI: 3.953757445491343
    },
    2021: {
      B11: 4.252402168270517,
      B12: 9.218157858862403,
      B4: 4.32725914239357,
      B8: 2.2749749285608,
      RBR: 9.935461889455542,
      RdNBR: 6.723836249659322,
      dNBR: 4.11103582700888,
      dNDVI: 2.810124408979414
    },
    2023: {
      B11: 8.677273150323858,
      B12: 11.715720457766512,
      B4: 10.734641332453869,
      B8: 13.725995816229037,
      RBR: 9.525462947212544,
      RdNBR: 11.460909237587034,
      dNBR: 9.999964941026638,
      dNDVI: 8.033658778211498
    }
  };

  // Ambil data importance berdasarkan tahun aktif
  var importance = importanceData[currentYear];
  if (!importance) {
    targetPanel.add(ui.Label('Data feature importance untuk tahun ' + currentYear + ' tidak tersedia.', {color: 'red'}));
    callback();
    return;
  }

  var fitur = Object.keys(importance);
  var nilai = fitur.map(function(k) { return [importance[k]] }); 

  // Buat dan tambahkan chart
  var chart = ui.Chart.array.values(nilai, 0, fitur)
  .setChartType('BarChart')
  .setOptions({
    title: 'Feature Importance Tahun ' + currentYear,
    height: 350, // Ukuran chart lebih besar
    width: 400,
    hAxis: {title: 'Fitur', textStyle: {fontSize: 12}},
    vAxis: {title: 'Nilai Importance', textStyle: {fontSize: 12}},
    series: {0: {color: '#1b9e77'}},
    legend: {position: 'none'},
    fontSize: 12,
    titleTextStyle: {fontSize: 14, bold: true}
  });

  targetPanel.add(chart);
  callback();
  }

  // --- Akurasi Data Validasi Lapangan ---
  function displayValidationAccuracy(targetPanel, callback) {
    targetPanel.add(ui.Label('Akurasi Validasi Lapangan:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
    var dataValidasi = {
      2019: 0.897959,
      2021: 0.938776,
      2023: 0.882653
    };
    var dataKappa = {
    2019: 0.608392,
    2021: 0.693111,
    2023: 0.401169
    };
     var confusionMatricesValidasi = {
    2019: ee.Array([[156,5], [15,20]]),
    2021: ee.Array([[168,4], [8,16]]),
    2023: ee.Array([[163,8], [15,10]])
    };
    var prodUserAccuracyValidasi = {
    2019: ee.Array([
      [0.9122, 0.9689],     // Kelas 0 (Tidak Terbakar)
      [0.8000, 0.571429]     // Kelas 1 (Terbakar)
    ]),
    2021: ee.Array([
      [0.9545, 0.9767],
      [0.8000, 0.666667]
    ]),
    2023: ee.Array([
      [0.9157, 0.9532],     
      [0.5555, 0.4000]
    ])
};
    // Ambil nilai berdasarkan currentYear
    var akurasi = dataValidasi[currentYear];
    var kappa = dataKappa[currentYear];
    var matrix = confusionMatricesValidasi[currentYear];
    var prodUserDict = prodUserAccuracyValidasi[currentYear];
    
    var formattedAkurasi = (akurasi !== null && akurasi !== undefined) ? akurasi.toFixed(4) : 'N/A';
    var formattedAkurasiPercent = (akurasi !== null && akurasi !== undefined) ? (akurasi * 100).toFixed(2) + '%' : 'N/A';
    var formattedKappa = (kappa !== null && kappa !== undefined) ? kappa.toFixed(4) : 'N/A';
    
    targetPanel.add(ui.Label(' • Overall Accuracy ' +  ': ' + formattedAkurasi+ ' (' + formattedAkurasiPercent + ')'));
    targetPanel.add(ui.Label(' • Kappa Coefficient  ' +  ': ' + formattedKappa ));
    targetPanel.add(ui.Label('Confusion Matrix:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
    
    function createConfMatrixChartValidasi(matrix, currentYear) {
    return ui.Chart.array.values({
        array: matrix,
        axis: 0,
        xLabels: ['Aktual: Tidak Terbakar (0)', 'Aktual: Terbakar (1)']
      })
      .setSeriesNames(['Prediksi: Tidak Terbakar (0)', 'Prediksi: Terbakar (1)'])
      .setChartType('Table')
      .setOptions({
        title: 'Confusion Matrix (Training ' + currentYear + ')',
        hAxis: {title: 'Prediksi'},
        vAxis: {title: 'Aktual'},
        colors: ['#cccccc', '#cccccc'],
        legend: {position: 'none'},
        fontSize: 10,
        height: 200,
        width: 350
      });
  }
  
  // Tambahkan confusion matrix hanya untuk tahun yang dipilih
  targetPanel.add(createConfMatrixChartValidasi(matrix, currentYear));

  targetPanel.add(ui.Label('Akurasi Spesifik:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
    
function createProdUserChartVal(prodUserDict, currentYear) {
  // Ambil array dan ubah ke persen
  var accuracyArray = prodUserDict[currentYear];
  var accuracyPercent = accuracyArray.multiply(100);

  // Buat label untuk menjelaskan isi tabel
  var explanation = ui.Label({
    value: 'Kolom 0: Producer Accuracy (%)\nKolom 1: User Accuracy (%)',
    style: { fontSize: '11px', color: 'gray', margin: '0 0 4px 0' }
  });

  // Buat chart tipe tabel
  var chart = ui.Chart.array.values({
      array: accuracyPercent,
      axis: 0,
      xLabels: ['Tidak Terbakar (0)', 'Terbakar (1)']
    })
    .setSeriesNames(['Akurasi Producer', 'Akurasi User'])
    .setChartType('Table')
    .setOptions({
      title: 'Producer & User Accuracy (Training) - ' + currentYear,
      fontSize: 12,
      legend: { position: 'none' },
      width: 360,
      height: 200
    });

  // Gabungkan label dan chart dalam satu panel
  return ui.Panel([explanation, chart]);
  }

// Contoh pemanggilan (pastikan 'produs' dan 'currentYear' sudah didefinisikan):
  targetPanel.add(createProdUserChartVal(prodUserAccuracyValidasi, currentYear));
    callback();
  }

  // --- Luas Terbakar ---
  function displayBurnedArea(targetPanel, callback) {
    targetPanel.add(ui.Label('Luas Terbakar:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
    // Data luas tetap per tahun (hektar)
  var dataLuasTetap = {
    2019: 47.68,
    2021: 23.64,
    2023: 45.68
  };

  // Ambil nilai berdasarkan currentYear
  var luas = dataLuasTetap[currentYear];
  var formattedLuas = (luas !== null && luas !== undefined) ? luas.toFixed(2) : '0.00';

  // Tampilkan ke panel
  targetPanel.add(ui.Label(' • Total Luas Terbakar Tahun ' + currentYear + ': ' + formattedLuas + ' hektar'));
  var linkLabel = ui.Label({
    value: '[Klik Disini: Tautan Tabel Hasil Analisis Luas Area Kebakaran Hutan dan Lahan ]',
    style: {fontSize: '12px', color: 'blue', textAlign: 'left', margin: '2px 0 2px 6px'},
    targetUrl: 'https://drive.google.com/file/d/1K-hlx8zEQE3pw2GTn7b9jYDQYwgpTsgm/view?usp=sharing'
  });
  targetPanel.add(linkLabel);
      callback();
  }

  // --- Grafik Perbandingan Luas Terbakar per Tahun ---
  function displayComparisonChart(targetPanel, callback) {
    targetPanel.add(ui.Label('📊 Titik Panas vs Luas Validasi Area Terbakar:',{fontWeight: 'bold', margin: '10px 0 5px 0'}));
    tampilkanGrafikTitikPanas(targetPanel);
    targetPanel.add(ui.Label('📊 Perbandingan Luas Terbakar per Tahun:', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
    tampilkanGrafikLuasTerbakar(targetPanel);// Panggil fungsi yang menambahkan ke targetPanel
    callback(); // Callback segera karena tampilkanGrafikLuasTerbakar sudah memiliki evaluasi internal
  }

  // =============================================================================
  // Tambahkan Tombol dan Panel Kontainer untuk Setiap Bagian (dalam urutan)
  // =============================================================================
  createToggleableSection(' Grafik ', displayComparisonChart);
  createToggleableSection('Luas Terbakar', displayBurnedArea);
  createToggleableSection(' Feature Importance', displayFeatureImportance);
  createToggleableSection(' Akurasi Training Model', displayTrainAccuracy);
  createToggleableSection(' Akurasi Testing Model', displayTestAccuracy);
  createToggleableSection(' Akurasi Validasi Lapangan', displayValidationAccuracy);

}