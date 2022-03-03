var resultsList = document.getElementById('results');
var response = null;
var query;
$( "#button" ).click(function() {
resultsList.innerHTML = "";
query = document.getElementById('search').value;
response = null;
load(query);
console.log("searching...");
});

function search(){
    var term = document.getElementById('search').value;
    queryJSON = {
        "term": term
    };
    $.ajax({
        url: '/search',
        type: 'POST',
        contentType: 'application/json; charset=utf-8',
        datatype: "json",
        data: JSON.stringify(queryJSON),
        success: function (data) {
            console.log("Recieving search results: " + JSON.stringify(data));
            response = data;
            renderResults(response);
            //document.getElementById("loading").remove();
        },
        error: function () {
            console.log('Error');
        }
    });
}

const timer = ms => new Promise(res => setTimeout(res, ms))
async function load(query){
    //tomorrow basically take the input of the text bar and use this to randomly generate a bunch of results after a brief loading animation and display them
    search();
    var loadingContent = document.createElement('div');
    loadingContent.innerHTML = '<p class="loading" id="loading">Loading.</p>'
    resultsList.prepend(loadingContent);
    var inc = 0;
    while(response == null){
        setTimeout(function(){
        if (response == null){
            loadingContent.innerHTML = '<p class="loading" id="loading">Loading..</p>';
        } else {
            document.getElementById("loading").remove();
        }
        }, 1000);
        setTimeout(function(){
        if (response == null){
            loadingContent.innerHTML = '<p class="loading" id="loading">Loading...</p>';
        } else {
            document.getElementById("loading").remove();
        }
        }, 2000);
        setTimeout(function(){
        if (response == null){
            loadingContent.innerHTML = '<p class="loading" id="loading">Loading</p>';
        } else {
            document.getElementById("loading").remove();
        }
        }, 3000);
        setTimeout(function(){
        if (response == null){
            loadingContent.innerHTML = '<p class="loading" id="loading">Loading.</p>';
        } else {
            document.getElementById("loading").remove();
        }
        }, 4000);
        await timer(4000);
        if (response != null){
            document.getElementById("loading").remove();
        }
    }
    //generateDemo(query);
}


function displayToggle(cmpdName,toggle){
    var props = document.getElementById(cmpdName+"-props");
    var biosig = document.getElementById(cmpdName+"-biosig");
    var assay = document.getElementById(cmpdName+"-assays");
    var synth = document.getElementById(cmpdName+"-synth");
    if(toggle == "props"){
        props.style.display = "block";
        biosig.style.display = "none";
        assay.style.display = "none";
        synth.style.display = "none";
    }
    if(toggle == "biosig"){
        props.style.display = "none";
        biosig.style.display = "block";
        assay.style.display = "none";
        synth.style.display = "none";
    }
    if(toggle == "assay"){
        props.style.display = "none";
        biosig.style.display = "none";
        assay.style.display = "block";
        synth.style.display = "none";
    }
    if(toggle == "synth"){
        props.style.display = "none";
        biosig.style.display = "none";
        assay.style.display = "none";
        synth.style.display = "block";
    }
}

function generateData(cols,bounds){
        let vals = [];
        for(let i = 0; i < cols; i++){
            vals.push(Math.random() * (bounds.max - bounds.min) + bounds.min);
        }
        return vals;
    }

function renderBiosig(cmpdName,inchi){
    var options = {
      series: [
      {
        name: 'Naloxone',
        data: generateData(9, {
          min: -30,
          max: 55
        })
      },
      {
        name: 'Desired biosignature',
        data: generateData(9, {
          min: -30,
          max: 55
        })
      }
    ],
      chart: {
      height: 350,
      type: 'heatmap',
    },
    plotOptions: {
      heatmap: {
        shadeIntensity: 0.5,
        radius: 0,
        useFillColorAsStroke: true,
        colorScale: {
          ranges: [{
              from: -30,
              to: 5,
              name: 'Partial agonism',
              color: '#00A100'
            },
            {
              from: 6,
              to: 20,
              name: 'Agonism',
              color: '#128FD9'
            },
            {
              from: 21,
              to: 45,
              name: 'Antagonism',
              color: '#FFB200'
            },
            {
              from: 46,
              to: 55,
              name: 'No effect',
              color: '#050505'
            }
          ]
        }
      }
    },
    dataLabels: {
      enabled: false
    },
    stroke: {
      width: 1
    },
    title: {
      text: cmpdName + ' — Last updated 1/1/2022'
    },
     xaxis: {
         type: 'category',
         tickPlacement: 'between',
         axisTicks: {
            show: true
         },
  tickAmount: 9,
  range: 9,
  labels: {
    show: true,
    rotate: -90,
    rotateAlways: true,
    formatter: function(value,timestamp,opts){
        console.log(value);
        return ["MOR","DOR","KOR","NOR","D2DD","D2LDR","D3DR","NMDAR","AMPAR"][value-1]
    }
  }
}
    };

    var chart = new ApexCharts(document.querySelector("#"+cmpdName+"-biosig"), options);
    chart.render();
}

function renderResults(results){
    var numEntries = results.ids.length;
    let entries = results.ids;
    let props = results.props;
    var res = document.getElementById('results');
    for(let i = 0; i < numEntries; i++){
        // LATER DO NOT INDEX INTO PROPS LIST; KEY IN BY INCHI KEY
        var cmpdName = entries[i].name;
        var inchikey = entries[i].inchikey;
        var smiles = entries[i].smiles;
        var MolecularFormula = props[i].pubchem.MolecularFormula;
        var inchi = entries[i].inchi;
        var molwt = props[i].pubchem.ExactMass;
        var hbonddonors = props[i].pubchem.HBondDonorCount;
        var hbondacceptors = props[i].pubchem.HBondAcceptorCount;
        var logp = props[i].rdkit.MolLogP;
        //var mol = RDKitModule.get_mol(smiles);
        //var svg = mol.get_svg();
        var svg = "";
        var entry = document.createElement("div");
        entry.innerHTML = '<div class="row" id="'+cmpdName+'" style="border-radius: 25px;border-color: var(--bs-indigo);"><div class="col-3 d-xl-flex justify-content-xl-end"><img src="./aspire_files/Screen Shot 2022-01-01 at 3.47.10 PM.png" width="200px" height="200px"></div><div class="col-9"><a style="color: grey;font-size: 30px;text-align: left;margin-bottom: 5px;" href="/summary/'+inchikey+'">'+cmpdName+'</a><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">SMILES<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+smiles+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">InChI Key<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+inchikey+'<br></p><div class="container"><div class="row"><div class="col-md-3"><button class="btn btn-primary border-white collapsed" type="button" style="width: 100%;background: var(--bs-gray-800);border-radius: 0px;border-style: none;" data-bs-toggle="collapse" data-bs-target="#'+cmpdName+'-collapse" aria-expanded="false" aria-controls="'+cmpdName+'-collapse" onclick="displayToggle(&#39;'+cmpdName+'&#39;,&#39;props&#39;)">Properties</button></div><div class="col-md-3"><button class="btn btn-primary border-white collapsed" type="button" style="width: 100%;background: var(--bs-gray-800);border-radius: 0px;border-style: none;" data-bs-toggle="collapse" data-bs-target="#'+cmpdName+'-collapse" aria-expanded="false" aria-controls="'+cmpdName+'-collapse" onclick="displayToggle(&#39;'+cmpdName+'&#39;,&#39;biosig&#39;)">Biosignature</button></div><div class="col-md-3"><button class="btn btn-primary border-white collapsed" type="button" style="width: 100%;background: var(--bs-gray-800);border-radius: 0px;border-style: none;" data-bs-toggle="collapse" data-bs-target="#'+cmpdName+'-collapse" aria-expanded="false" aria-controls="'+cmpdName+'-collapse" onclick="displayToggle(&#39;'+cmpdName+'&#39;,&#39;assay&#39;)">Assays</button></div><div class="col-md-3"><button class="btn btn-primary border-white collapsed" type="button" style="width: 100%;background: var(--bs-gray-800);border-radius: 0px;border-style: none;" data-bs-toggle="collapse" data-bs-target="#'+cmpdName+'-collapse" aria-expanded="false" aria-controls="'+cmpdName+'-collapse" onclick="displayToggle(&#39;'+cmpdName+'&#39;,&#39;synth&#39;)">Synthesis</button></div></div></div><div class="row" style="border-style: none;"><div class="col" style="border-style: none;"><div class="row"><div class="col"><div id="'+cmpdName+'-collapse" class="collapse" style="padding: 10px;"><div id="'+cmpdName+'-biosig" style="min-height: 365px; display: none;"><div id="apexchartsqfpbc9xf" class="apexcharts-canvas apexchartsqfpbc9xf" style="width: 0px; height: 350px;"><svg id="SvgjsSvg1181" width="0" height="350" xmlns="http://www.w3.org/2000/svg" version="1.1" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:svgjs="http://svgjs.dev" class="apexcharts-svg" xmlns:data="ApexChartsNS" transform="translate(0, 0)" style="background: transparent;"><g id="SvgjsG1184" class="apexcharts-annotations"></g><g id="SvgjsG1183" class="apexcharts-inner apexcharts-graphical"><defs id="SvgjsDefs1182"></defs></g></svg><div class="apexcharts-legend"></div></div></div><div id="'+cmpdName+'-props" style="display: block;"><div class="row"><div class="col"><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">InChI<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+inchi+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Molecular Formula<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+MolecularFormula+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">CAS<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+'Not returned'+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Molecular Weight<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+molwt+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Hydrogen Bond Donor Count<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+hbonddonors+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Hydrogen Bond Acceptor Count<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+hbondacceptors+'<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">logP<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">'+logp+'<br></p></div></div></div><div id="'+cmpdName+'-assays"><div class="row"><div class="col"></div></div></div><div id="'+cmpdName+'-synth"></div></div></div></div></div></div></div><div class="col"><hr></div></div>'
        res.appendChild(entry);
        
        renderBiosig(cmpdName,"");
    }
    
    return vals;
}


function generateDemo(query){
    var numEntries = Math.floor(Math.random() * (10 - 5) + 5)
    let entries = [];
    var res = document.getElementById('results');
    for(let i = 0; i < numEntries; i++){
        var cmpdName = query + String(Math.floor(Math.random() * (10 - 5) + 5));
        var entry = document.createElement("div");
        entry.innerHTML = '<div class="row" id="'+cmpdName+'" style="border-radius: 25px;border-color: var(--bs-indigo);"><div class="col-3 d-xl-flex justify-content-xl-end"><img src="./aspire_files/Screen Shot 2022-01-01 at 3.47.10 PM.png" width="200px" height="200px"></div><div class="col-9"><p style="color: grey;font-size: 30px;text-align: left;margin-bottom: 5px;">'+cmpdName+'</p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">SMILES<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">C=CCN1CCC23C4C(=O)CCC2(C1CC5=C3C(=C(C=C5)O)O4)O<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">InChI Key<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">UZHSEJADLWPNLE-GRGSLBFTSA-N<br></p><div class="container"><div class="row"><div class="col-md-3"><button class="btn btn-primary border-white collapsed" type="button" style="width: 100%;background: var(--bs-gray-800);border-radius: 0px;border-style: none;" data-bs-toggle="collapse" data-bs-target="#'+cmpdName+'-collapse" aria-expanded="false" aria-controls="'+cmpdName+'-collapse" onclick="displayToggle(&#39;'+cmpdName+'&#39;,&#39;props&#39;)">Properties</button></div><div class="col-md-3"><button class="btn btn-primary border-white collapsed" type="button" style="width: 100%;background: var(--bs-gray-800);border-radius: 0px;border-style: none;" data-bs-toggle="collapse" data-bs-target="#'+cmpdName+'-collapse" aria-expanded="false" aria-controls="'+cmpdName+'-collapse" onclick="displayToggle(&#39;'+cmpdName+'&#39;,&#39;biosig&#39;)">Biosignature</button></div><div class="col-md-3"><button class="btn btn-primary border-white collapsed" type="button" style="width: 100%;background: var(--bs-gray-800);border-radius: 0px;border-style: none;" data-bs-toggle="collapse" data-bs-target="#'+cmpdName+'-collapse" aria-expanded="false" aria-controls="'+cmpdName+'-collapse" onclick="displayToggle(&#39;'+cmpdName+'&#39;,&#39;assay&#39;)">Assays</button></div><div class="col-md-3"><button class="btn btn-primary border-white collapsed" type="button" style="width: 100%;background: var(--bs-gray-800);border-radius: 0px;border-style: none;" data-bs-toggle="collapse" data-bs-target="#'+cmpdName+'-collapse" aria-expanded="false" aria-controls="'+cmpdName+'-collapse" onclick="displayToggle(&#39;'+cmpdName+'&#39;,&#39;synth&#39;)">Synthesis</button></div></div></div><div class="row" style="border-style: none;"><div class="col" style="border-style: none;"><div class="row"><div class="col"><div id="'+cmpdName+'-collapse" class="collapse" style="padding: 10px;"><div id="'+cmpdName+'-biosig" style="min-height: 365px; display: none;"><div id="apexchartsqfpbc9xf" class="apexcharts-canvas apexchartsqfpbc9xf" style="width: 0px; height: 350px;"><svg id="SvgjsSvg1181" width="0" height="350" xmlns="http://www.w3.org/2000/svg" version="1.1" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:svgjs="http://svgjs.dev" class="apexcharts-svg" xmlns:data="ApexChartsNS" transform="translate(0, 0)" style="background: transparent;"><g id="SvgjsG1184" class="apexcharts-annotations"></g><g id="SvgjsG1183" class="apexcharts-inner apexcharts-graphical"><defs id="SvgjsDefs1182"></defs></g></svg><div class="apexcharts-legend"></div></div></div><div id="'+cmpdName+'-props" style="display: block;"><div class="row"><div class="col"><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">InChI<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">InChI=1S/C19H21NO4/c1-2-8-20-9-7-18-15-11-3-4-12(21)16(15)24-17(18)13(22)5-6-19(18,23)14(20)10-11/h2-4,14,17,21,23H,1,5-10H2/t14-,17+,18+,19-/m1/s1<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Molecular Formula<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">CAS<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Molecular Weight<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Hydrogen Bond Donor Count<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Hydrogen Bond Acceptor Count<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">logP</p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Bioassay Link</p></div></div></div><div id="'+cmpdName+'-assays"><div class="row"><div class="col"><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">InChI<br></p><p style="color: grey;font-size: 14px;text-align: left;margin-bottom: 10px;">InChI=1S/C19H21NO4/c1-2-8-20-9-7-18-15-11-3-4-12(21)16(15)24-17(18)13(22)5-6-19(18,23)14(20)10-11/h2-4,14,17,21,23H,1,5-10H2/t14-,17+,18+,19-/m1/s1<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Molecular Formula<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">CAS<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Molecular Weight<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Hydrogen Bond Donor Count<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Hydrogen Bond Acceptor Count<br></p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">logP</p><p style="color: grey;font-size: 20px;text-align: left;margin-bottom: 0px;">Bioassay Link</p></div></div></div><div id="'+cmpdName+'-synth"></div></div></div></div></div></div></div><div class="col"><hr></div></div>'
        res.appendChild(entry);
        renderBiosig(cmpdName,"");
    }
    
    return vals;
}
