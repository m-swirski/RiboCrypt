(function(){
  const colors={'Mammals':'#c9563c','Other vertebrates':'#e69b4d','Invertebrates':'#d4b549','Plants':'#55a66c','Fungi':'#9476b7','Protists':'#5aada6','Bacteria':'#5388c1','Archaea':'#98a8bd','Viruses':'#b7b9b1','Unclassified':'#c8d0cb'};
  function draw(d){
    if(!Array.isArray(d)) d=(d.name||[]).map((name,i)=>({name,clade:d.clade[i],libraries:d.libraries[i],studies:d.studies[i]}));
    // The wheel depicts species with reported libraries. Keep zero-count
    // organisms in the separate species selector, where they remain usable.
    d=d.filter(s=>Number(s.libraries)>0);
    const svg=document.getElementById('rc-home-wheel'),legend=document.getElementById('rc-home-legend'),tip=document.getElementById('rc-home-tooltip');if(!svg)return;
    svg.replaceChildren();legend.replaceChildren();const ns='http://www.w3.org/2000/svg';
    const node=(n,a)=>{let e=document.createElementNS(ns,n);Object.entries(a).forEach(([k,v])=>e.setAttribute(k,v));svg.appendChild(e);return e};
    const n=d.length,max=Math.max(1,...d.map(s=>s.libraries||0));
    const groups=[...new Set(d.map(s=>s.clade))],gap=1.8,slots=n+groups.length*gap;
    node('circle',{cx:310,cy:310,r:142,fill:'#f8fbf9',stroke:'#eaf0ed'});
    const point=(r,a)=>[310+r*Math.cos(a),310+r*Math.sin(a)];
    d.forEach((s,i)=>{const a=-Math.PI/2+2*Math.PI*(i+(groups.indexOf(s.clade)+.5)*gap)/slots,b=a+2*Math.PI/slots*.78,r=175+95*Math.log1p(s.libraries||0)/Math.log1p(max),A=point(153,a),B=point(r,a),C=point(r,b),D=point(153,b);
      const p=node('path',{d:`M ${A} L ${B} A ${r} ${r} 0 0 1 ${C} L ${D} A 153 153 0 0 0 ${A} Z`,fill:colors[s.clade]||colors.Unclassified,tabindex:0,role:'button','aria-label':`${s.name}: ${s.libraries} libraries, ${s.studies} studies`});
      const choose=()=>Shiny.setInputValue('home_wheel_species',s.name,{priority:'event'});p.addEventListener('click',choose);p.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();choose()}});
      const show=()=>{tip.textContent=`${s.name} · ${s.clade} · ${s.libraries} libraries · ${s.studies} studies`;tip.style.display='block';const box=p.getBoundingClientRect(),parent=tip.parentNode.getBoundingClientRect();tip.style.left=(box.left-parent.left)+'px';tip.style.top=(box.top-parent.top-35)+'px'};p.addEventListener('mouseenter',show);p.addEventListener('focus',show);p.addEventListener('mouseleave',()=>tip.style.display='none');p.addEventListener('blur',()=>tip.style.display='none');
    });
    [['EXPLORE THE ATLAS',286,11],[(n||0)+' species',320,23],[(d.reduce((a,s)=>a+(s.libraries||0),0)).toLocaleString()+' Ribo-seq libraries',350,11]].forEach(([t,y,size])=>{let x=node('text',{x:310,y,'text-anchor':'middle','font-size':size,fill:'#366058','pointer-events':'none'});x.textContent=t});
    Object.entries(colors).forEach(([group,color])=>{if(!d.some(x=>x.clade===group))return;let s=document.createElement('span'),i=document.createElement('i');i.style.background=color;s.append(i,document.createTextNode(group));legend.append(s)});
  }
  Shiny.addCustomMessageHandler('rc-home-wheel',draw);
})();
