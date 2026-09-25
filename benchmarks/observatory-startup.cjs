// Fresh sessions opening Select libraries. PROFILE=1 instruments Plotly separately.
const {chromium}=require('playwright');
const fs=require('node:fs');
const assert=require('node:assert/strict');
const urls=(process.argv[2] || 'http://127.0.0.1:4723/#Observatory').split(',');
const visits=Number(process.argv[3] || 6);

function instrument(profile) {
  const record=window.obsTiming={values:[],draws:[],calls:[]};
  if(profile) {
    let plotly;
    Object.defineProperty(window,'Plotly',{configurable:true,get:()=>plotly,set:value=>{
      plotly=value;
      for(const name of ['newPlot','react','restyle','relayout','addTraces','deleteTraces']) {
        const original=value[name];
        value[name]=function(...args) {
          const el=typeof args[0]==='string' ? document.getElementById(args[0]):args[0];
          const call={name,id:el?.id,start:performance.now(),
            patch:name==='relayout'?args[1]:undefined,
            caller:new Error().stack.split('\n').slice(2,5).join('\n')}; record.calls.push(call);
          const result=original.apply(this,args);call.syncEnd=performance.now();
          if(result?.then) result.then(()=>{call.end=performance.now();});
          return result;
        };
      }
    }});
  }
  const timer=setInterval(()=>{
    if(!window.jQuery) return;
    clearInterval(timer);
    jQuery(document).on('shiny:connected',()=>{record.connected=performance.now();});
    jQuery(document).on('shiny:value',event=>{
      record.values.push({id:event.name,time:performance.now()});
    });
    jQuery(document).on('draw.dt',event=>{
      if(event.target.closest('#observatory-selector-libraries_data_table')) record.draws.push(performance.now());
    });
  },1);
}

(async()=>{
  const browser=await chromium.launch({executablePath:process.env.CHROME || '/opt/google/chrome/chrome',
    headless:true,args:['--no-sandbox']});
  const results=[];
  try {
    for(let trial=0;trial<visits*urls.length;trial++) {
      const visit=Math.floor(trial/urls.length), url=urls[trial%urls.length];
      const context=await browser.newContext({viewport:{width:1440,height:1000}});
      const page=await context.newPage(), errors=[];
      page.on('pageerror',error=>errors.push(error.message));
      await page.addInitScript(instrument,Boolean(process.env.PROFILE));
      await page.goto(url,{waitUntil:'domcontentloaded'});
      await page.waitForFunction(()=>{
        const r=window.obsTiming, prefix='observatory-selector-';
        const umap=document.getElementById(prefix+'libraries_umap_plot');
        if(umap?._fullData?.some(t=>t.x?.length) && umap.getBoundingClientRect().width>0 &&
           umap.querySelector('.main-svg')) r.umap ||= performance.now();
        const table=document.querySelector('#'+prefix+'libraries_data_table .dataTables_scrollBody table');
        if(table && window.jQuery?.fn.dataTable?.isDataTable(table)) {
          const api=jQuery(table).DataTable(), info=api.page.info();
          if(info.recordsTotal>0 && api.rows({page:'current'}).data().length>0 &&
             table.getBoundingClientRect().width>0) r.dt ||= performance.now();
        }
        if(!r.umap || !r.dt) return false;
        r.both ||= performance.now();
        if(document.documentElement.classList.contains('shiny-busy')) return false;
        r.idle ||= performance.now();return true;
      },null,{timeout:120000,polling:25});
      await page.waitForTimeout(1500);
      const result=await page.evaluate(()=>{
        const prefix='observatory-selector-', el=document.getElementById(prefix+'libraries_umap_plot');
        const table=document.querySelector('#'+prefix+'libraries_data_table .dataTables_scrollBody table');
        return {...window.obsTiming,responseEnd:performance.getEntriesByType('navigation')[0].responseEnd,
          collection:document.getElementById(prefix+'dff')?.value,
          points:el.data.reduce((n,t)=>n+(t.x?.length || 0),0),traces:el.data.length,
          table:jQuery(table).DataTable().page.info(),
          rows:jQuery(table).DataTable().rows({page:'current'}).data().length,
          resources:performance.getEntriesByType('resource').map(r=>({name:r.name,start:r.startTime,
            end:r.responseEnd,bytes:r.transferSize,initiator:r.initiatorType})),
          errors:[...document.querySelectorAll('.shiny-output-error')].map(e=>e.textContent)};
      });
      assert.deepEqual(errors,[]);assert.deepEqual(result.errors,[]);
      assert.equal(result.collection,'all_samples-Homo_sapiens');
      assert.equal(result.rows,15);assert.ok(result.table.recordsTotal>1000);
      results.push({url,visit,...result});
      console.log(JSON.stringify({url,visit,...result,resources:undefined,calls:process.env.PROFILE?result.calls:undefined}));
      if(visit===0) await page.screenshot({path:'/tmp/observatory-startup-'+new URL(url).port+'.png',fullPage:true});
      await context.close();
    }
    fs.writeFileSync(process.env.RESULTS || '/tmp/observatory-startup.json',JSON.stringify(results,null,2));
  } finally {await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
