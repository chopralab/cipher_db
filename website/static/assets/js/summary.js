function generateData(cols,bounds){
    let vals = [];
    for(let i = 0; i < cols; i++){
        vals.push(Math.random() * (bounds.max - bounds.min) + bounds.min);
    }
    return vals;
}
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
      text: ''
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

    var chart = new ApexCharts(document.querySelector("#sig"), options);
    chart.render();
