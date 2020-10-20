
    var chart_config = {
        chart: {
            container: "#aepl",
            //rootOrientation: "WEST",
            connectors: {
              type: 'curve'
              //type: 'bCurve'
              //type: 'step'
              //type: 'straight'
            },
            animateOnInit: true,

            node: {
                collapsable: true
            },
            animation: {
                //nodeAnimation: "easeOutBounce",
                nodeSpeed: 700,
                //connectorsAnimation: "bounce",
                connectorsSpeed: 400
            }
        },

        nodeStructure: {
          text: { name: "aepl" },
          HTMLclass : 'dirNode',

          children: [




            {

              pseudo: true,
              //childrenDropLevel: 1,
              children: [
                {
                  text: { name: "aeplUtil" },
                  HTMLclass : 'dirNode',
                  //stackChildren: true,
                    //childrenDropLevel: 1,
                    children: [
                      {
                        text: { name: "buildcmap.m" },
                        HTMLclass : 'fileNode',
                        collapsed: true,
                      },
                      {
                        text: { name: "makeConditDict.m" },
                        HTMLclass : 'fileNode',
                        collapsed: true,
                      },
                    ]
                }
              ]

            },
            {
              pseudo: true,
              childrenDropLevel: 1,
              children: [
                {
                  text: { name: "processCzi" },
                  HTMLclass : 'dirNode',

                  children: [
                    {
                      text: { name: "(Run.m)", title: "(Run2.m)" },
                      HTMLclass : 'fileNode',
                      //collapsed: true,
                      children: [
                        {
                          text: { name: "(Mak2eConditDict.m)" },
                          HTMLclass : 'fileNode',
                          collapsed: true,
                          children: [
                            {
                              text: { desc: "idk"  },
                              HTMLclass : 'infoNode',
                            }
                          ]
                        },
                        {
                          text: { name: "<MicroscopeData.", title: "Original.", desc: "ReadData>" },
                          HTMLclass : 'fileNode',
                          collapsed: true,
                        },
                        {
                          text: { name: "segIms.m" },
                          HTMLclass : 'fileNode',
                          //collapsed: true,
                          children: [
                            {
                              text: { name: "segTexture_MSKCC.m" },
                              HTMLclass : 'fileNode',
                              collapsed: true,
                            },


                          ]
                        },
                        // {
                        //   pseudo: true,
                        //   //childrenDropLevel: ,
                        //   children: [
                            {
                              text: { name: "getTracks.m" },
                              HTMLclass : 'fileNode',
                              collapsed: true,
                              children: [
                                {
                                  text: { desc: "idk"  },
                                  HTMLclass : 'infoNode',
                                }
                              ]
                            },
                            {
                              text: { name: "DrawTracks.m", title: "ExportTrackStats.m" },
                              HTMLclass : 'fileNode',
                              //collapsed: true,
                              children: [
                                {
                                  text: { name: "MakeMovie.m" },
                                  HTMLclass : 'fileNode',
                                  collapsed: true,
                                  children: [
                                    {
                                      text: { desc: "idk"  },
                                      HTMLclass : 'infoNode',
                                    }
                                  ]
                                },
                              ]
                            },
                          // ]
                        // },


                      ]
                    },








                  ]
                }
              ]

          },
          {
            pseudo: true,
            childrenDropLevel: 4,
            children: [
                {
                  text: { name: "readPlotCsv" },
                  HTMLclass : 'dirNode',
                  //stackChildren: true,
                  children: [
                    {
                      text: { name: "plotCondits.m" },
                      HTMLclass : 'fileNode',
                      collapsed: true,
                    },
                    {
                      text: { name: "readd.m" },
                      HTMLclass : 'fileNode',
                      collapsed: true,
                    },
                    {
                      text: { name: "readPlateMap.m" },
                      HTMLclass : 'fileNode',
                      collapsed: true,
                    },
                    {
                      text: { name: "run.m" },
                      HTMLclass : 'fileNode',
                      collapsed: true,
                    },
                    {
                      text: { name: "scatDisp.m" },
                      HTMLclass : 'fileNode',
                      collapsed: true,
                    },

                    {
                      pseudo: true,
                      childrenDropLevel: 2,
                      children: [
                        {
                          text: { name: "scatDispAdd.m" },
                          HTMLclass : 'fileNode',
                          collapsed: true,
                        },
                        {
                          text: { name: "tempFigPdfSize.m" },
                          HTMLclass : 'fileNode',
                          collapsed: true,
                        },
                        {
                          text: { name: "Untitled.m" },
                          HTMLclass : 'fileNode',
                          collapsed: true,
                        },
                        {
                          text: { name: "wellConditAvgs.m" },
                          HTMLclass : 'fileNode',
                          collapsed: true,
                        },
                        {
                          text: { name: "wellConditMeds.m" },
                          HTMLclass : 'fileNode',
                          collapsed: true,
                        }
                      ]
                    }
                  ]



                }
              ]
            }


          ]
        }
    };
