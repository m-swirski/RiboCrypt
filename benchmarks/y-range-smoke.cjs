const {chromium}=require('playwright');
const assert=require('node:assert/strict');
const base=process.argv[2] || 'http://127.0.0.1:7821/';

async function generate(page,prefix) {
  await page.locator('#'+prefix+'toggle_settings').click();
  await page.locator('#'+prefix+'go').click();
}

async function range(page,id,expected) {
  await page.waitForFunction(({id,expected})=>{
    const plot=document.getElementById(id);
    return plot?._fullLayout && !document.documentElement.classList.contains('shiny-busy') &&
      JSON.stringify(plot._fullLayout.yaxis.range)===JSON.stringify(expected);
  },{id,expected},{timeout:120000});
}

(async()=>{
  const browser=await chromium.launch({executablePath:'/opt/google/chrome/chrome',
    headless:true,args:['--no-sandbox']});
  try {
    for(const obs of [false,true]) {
      const context=await browser.newContext({viewport:{width:1440,height:1000}});
      const page=await context.newPage(), errors=[];
      page.on('pageerror',e=>errors.push(e.message));
      const p=obs?'observatory-browser_obs-':'browser-', id=p+(obs?'browser_plot':'c');
      let copied;
      page.on('websocket',s=>s.on('framereceived',event=>{
        try {copied=JSON.parse(event.payload).custom?.['ribocrypt-copy-url']?.text || copied;} catch {}
      }));
      await page.goto(base+(obs?'#Observatory':'#browser'));
      if(obs) {
        await page.locator('a[data-value="Browse"]').click();
        await page.waitForFunction(p=>document.getElementById(p+'tx')?.value,p,{timeout:120000});
        await page.locator('#'+p+'go').click();
      }
      await page.waitForFunction(id=>document.getElementById(id)?._fullData?.length>0,id,{timeout:120000});
      const initial=await page.evaluate(id=>{
        const el=document.getElementById(id);
        return {range:el._fullLayout.yaxis.range,bottom:el._fullLayout.yaxis2.range};
      },id);
      await page.locator('#'+p+'toggle_settings').click();
      await page.locator('#'+p+'floating_settings a[data-value="Settings"]').click();
      assert.equal(await page.locator('#'+p+'y_range').inputValue(),'auto');
      await page.locator('#'+p+'y_range').fill('500k');
      await generate(page,p);
      await range(page,id,[0,500000]);
      await page.locator('#'+p+'toggle_settings').click();
      assert.deepEqual(await page.evaluate(id=>document.getElementById(id)._fullLayout.yaxis2.range,id),initial.bottom);
      await page.locator('#'+p+'clip_button').click();
      for(let i=0;!copied&&i<100;i++) await page.waitForTimeout(100);
      assert.ok(copied);
      const second=await context.newPage();second.on('pageerror',e=>errors.push(e.message));
      await second.goto(copied);await range(second,id,[0,500000]);
      await second.close();
      await page.locator('#'+p+'y_range').fill('100:20');
      await generate(page,p);
      await page.getByText('Invalid Y-axis range',{exact:true}).waitFor();
      assert.match(await page.locator('.modal-body').innerText(),/maximum greater than minimum/);
      await page.locator('.modal-footer button').click();
      await page.locator('#'+p+'toggle_settings').click();
      await page.locator('#'+p+'y_range').fill('10:600000');
      await generate(page,p);await range(page,id,[10,600000]);
      await page.locator('#'+p+'toggle_settings').click();
      await page.screenshot({path:'/tmp/y-range-'+(obs?'observatory':'browser')+'.png'});
      await page.locator('#'+p+'y_range').fill('auto');
      await generate(page,p);await range(page,id,initial.range);
      assert.deepEqual(errors,[]);
      assert.deepEqual(await page.locator('.shiny-output-error').allTextContents(),[]);
      console.log(JSON.stringify({obs,initial:initial.range,fixed:true,url:true,invalidRecovery:true,auto:true,errors}));
      await context.close();
    }
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
