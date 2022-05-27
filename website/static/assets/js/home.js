function getList(){
    queryJSON = {};
    $.ajax({
        url: '/top',
        type: 'POST',
        contentType: 'application/json; charset=utf-8',
        datatype: "json",
        data: JSON.stringify(queryJSON),
        success: function (data) {
            console.log("Recieving search results: " + JSON.stringify(data));
            response = data;
            console.log(response)
            renderResults(response);
            //document.getElementById("loading").remove();
        },
        error: function () {
            console.log('Error');
        }
    });
}

getList();

function generateData(cols,bounds){
    let vals = [];
    for(let i = 0; i < cols; i++){
        vals.push(Math.random() * (bounds.max - bounds.min) + bounds.min);
    }
    return vals;
}

function renderBiosigs(desired, top){
    let desList = [];
    let t1 = [];
    let t2 = [];
    let t3 = [];
    let t4 = [];
    let t5 = [];
    desList.push(desired[1].AMPAR);
    desList.push(desired[1].D2LDR);
    desList.push(desired[1].DRD2);
    desList.push(desired[1].DRD3);
    desList.push(desired[1].NMDAR);
    desList.push(desired[1].deltaOR);
    desList.push(desired[1].kappaOR);
    desList.push(desired[1].muOR);
    desList.push(desired[1].nociceptinOR);
    t1.push(top[0][1].AMPAR);
    t1.push(top[0][1].D2LDR);
    t1.push(top[0][1].DRD2);
    t1.push(top[0][1].DRD3);
    t1.push(top[0][1].NMDAR);
    t1.push(top[0][1].deltaOR);
    t1.push(top[0][1].kappaOR);
    t1.push(top[0][1].muOR);
    t1.push(top[0][1].nociceptinOR);
    t2.push(top[1][1].AMPAR);
    t2.push(top[1][1].D2LDR);
    t2.push(top[1][1].DRD2);
    t2.push(top[1][1].DRD3);
    t2.push(top[1][1].NMDAR);
    t2.push(top[1][1].deltaOR);
    t2.push(top[1][1].kappaOR);
    t2.push(top[1][1].muOR);
    t2.push(top[1][1].nociceptinOR);
    t3.push(top[2][1].AMPAR);
    t3.push(top[2][1].D2LDR);
    t3.push(top[2][1].DRD2);
    t3.push(top[2][1].DRD3);
    t3.push(top[2][1].NMDAR);
    t3.push(top[2][1].deltaOR);
    t3.push(top[2][1].kappaOR);
    t3.push(top[2][1].muOR);
    t3.push(top[2][1].nociceptinOR);
    t4.push(top[3][1].AMPAR);
    t4.push(top[3][1].D2LDR);
    t4.push(top[3][1].DRD2);
    t4.push(top[3][1].DRD3);
    t4.push(top[3][1].NMDAR);
    t4.push(top[3][1].deltaOR);
    t4.push(top[3][1].kappaOR);
    t4.push(top[3][1].muOR);
    t4.push(top[3][1].nociceptinOR);
    t5.push(top[4][1].AMPAR);
    t5.push(top[4][1].D2LDR);
    t5.push(top[4][1].DRD2);
    t5.push(top[4][1].DRD3);
    t5.push(top[4][1].NMDAR);
    t5.push(top[4][1].deltaOR);
    t5.push(top[4][1].kappaOR);
    t5.push(top[4][1].muOR);
    t5.push(top[4][1].nociceptinOR);

    
    var options = {
          series: [{
            name: 'NAA Candidate 2575',
            data: t5
          },
          {
            name: 'NAA Candidate 3180',
            data: t4
          },
          {
            name: 'NAA Candidate 4442',
            data: t3
          },
          {
            name: 'NAA Candidate 4441',
            data: t2
          },
          {
            name: 'NAA Candidate 516',
            data: t1
          },
          {
            name: 'NAA Control Signature',
            data: desList
          }
        ],
          chart: {
          height: 350,
          type: 'heatmap',
        },
        plotOptions: {
          heatmap: {
            shadeIntensity: 1,
            radius: 0,
            useFillColorAsStroke: true,
            colorScale: {
              ranges: [{
                  from: 0,
                  to: 1,
                  name: 'Partial agonism',
                  color: '#00A100'
                },
                {
                  from: 1,
                  to: 2,
                  name: 'Agonism',
                  color: '#128FD9'
                },
                {
                  from: 2,
                  to: 3,
                  name: 'Antagonism',
                  color: '#FFB200'
                },
                {
                  from: 3,
                  to: 4,
                  name: 'Unknown effect',
                  color: '#050505'
                },
                {
                  from: -9999999999999999999,
                  to: -1,
                  name: 'Null',
                  color: '#ed3e3e'
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
          text: 'Top hits'
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
            return ["AMPAR","D2LDR","DRD2","DRD3","NMDAR","deltaOR","kappaOR","muOR","nociceptinOR"][value-1]
        }
      }
    }
        };

        var chart = new ApexCharts(document.querySelector("#biosig"), options);
        chart.render();
}




function renderResults(results){
    console.log(results.top);
    console.log(results.desired);
    let desired = results.desired;
    let top = results.top;
    for (var key of Object.keys(desired[1])){
        if (desired[1][key] == -1){
            continue;
        }
        if (desired[2][key] == "partial agonist"){
            desired[1][key] = desired[1][key];
        }
        if (desired[2][key] == "binder"){
            desired[1][key] = desired[1][key] + 3;
        }
        if (desired[2][key] == "agonism"){
            desired[1][key] = desired[1][key] + 1;
        }
        if (desired[2][key] == "antagonism"){
            desired[1][key] = desired[1][key] + 2;
        }
        if (desired[2][key] == "unknown"){
            desired[1][key] = desired[1][key] + 3;
        }
    }
    
    for (var i=0; i<5; i++){
        for (var key of Object.keys(top[i][1])){
            if (top[i][1][key] == -1){
                continue;
            }
            if (top[i][2][key] == "partial agonist"){
                top[i][1][key] = top[i][1][key];
            }
            if (top[i][2][key] == "binder"){
                top[i][1][key] = top[i][1][key] + 3;
            }
            if (top[i][2][key] == "agonism"){
                top[i][1][key] = top[i][1][key] + 1;
            }
            if (top[i][2][key] == "antagonism"){
                top[i][1][key] = top[i][1][key] + 2;
            }
            if (top[i][2][key] == "unknown"){
                top[i][1][key] = top[i][1][key] + 3;
            }
        }
    }
    console.log(results.top);
    console.log(results.desired);
    renderBiosigs(desired,top);
}
