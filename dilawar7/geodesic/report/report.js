/* SPDX-License-Identifier: Apache-2.0 */
(() => {
  'use strict';
  const data = JSON.parse(document.getElementById('path-data').textContent);
  const $ = id => document.getElementById(id);
  const count = data.metadata.nimages;
  const state = {frame:0, source:'geodesic', timer:null, active:'dihedral_1_2_3_4', picking:false, picks:[]};
  const selected = {distance:'distance_1_4', angle:'angle_1_2_3', dihedral:'dihedral_1_2_3_4'};
  const names = {distance:'Distance', angle:'Angle', dihedral:'Dihedral'};
  let viewer = null;
  let initialView = null;
  const chartParts = {};
  const entryById = id => data.measurements.find(entry => entry.id === id);
  const atomLabel = atom => data.symbols[atom] + (atom + 1);
  const textValue = (value, unit) => value === null || !Number.isFinite(value) ? 'Undefined' : value.toFixed(unit === 'Å' ? 3 : 1) + ' ' + unit;

  function renderMolecule(reset=false) {
    if (!viewer) return;
    const view = viewer.getView();
    viewer.removeAllModels(); viewer.removeAllLabels(); viewer.removeAllShapes();
    const xyz = data.coordinates[state.source][state.frame];
    const active = entryById(state.active);
    const highlighted = state.picking ? state.picks : active.atoms;
    const atoms = xyz.map((point, i) => {
      const adjacent = data.bonds.filter(pair => pair.includes(i)).map(pair => pair[0] === i ? pair[1] : pair[0]);
      return {elem:data.symbols[i], x:point[0], y:point[1], z:point[2],
              serial:i+1, index:i, bonds:adjacent, bondOrder:adjacent.map(() => 1)};
    });
    const model = viewer.addModel(); model.addAtoms(atoms);
    const carbon = {sphere:{scale:.26,color:'#475569'},stick:{radius:.14,color:'#7b8a9d'}};
    const hydrogen = {sphere:{scale:.22,color:'#d9e2ed'},stick:{radius:.07,color:'#b4c1d0'}};
    model.setStyle({elem:'C'}, carbon);
    model.setStyle({elem:'H'}, $('show-h').checked ? hydrogen : {sphere:{hidden:true},stick:{hidden:true}});
    model.setStyle({index:highlighted.filter(i => data.symbols[i] !== 'H' || $('show-h').checked)},
                   {sphere:{scale:.29,color:'#087f8c'},stick:{radius:.15,color:'#087f8c'}});
    if ($('show-labels').checked) {
      atoms.forEach((atom,i) => {
        if (data.symbols[i] === 'C' || (highlighted.includes(i) && $('show-h').checked)) {
          viewer.addLabel(atomLabel(i), {position:atom, fontSize:14, fontColor:0x172944,
                          fontOpacity:1, showBackground:true, backgroundColor:0xffffff,
                          backgroundOpacity:.85, padding:2, inFront:true, screenOffset:{x:12,y:12}});
        }
      });
    }
    if (highlighted.length >= 2) {
      for (let i=1;i<highlighted.length;i++) {
        const point = index => ({x:xyz[index][0],y:xyz[index][1],z:xyz[index][2]});
        viewer.addLine({start:point(highlighted[i-1]),end:point(highlighted[i]),
                        color:'#087f8c',linewidth:2,dashed:true});
      }
    }
    model.setClickable({}, true, atom => {
      const index = atom.index;
      if (!state.picking) {
        $('custom-status').textContent = atomLabel(index) + ' · atom number ' + (index+1) + '. Use Pick atoms to define a measurement.';
        return;
      }
      if (!state.picks.includes(index)) state.picks.push(index);
      $('custom-atoms').value = state.picks.map(i => i+1).join(', ');
      const required = {distance:2,angle:3,dihedral:4}[$('custom-kind').value];
      if (state.picks.length === required) $('custom-form').requestSubmit();
      else renderMolecule();
    });
    if (!initialView || reset) {
      viewer.zoomTo(); viewer.rotate(18,'y'); viewer.rotate(15,'x');
      initialView = viewer.getView();
    } else viewer.setView(view);
    viewer.render();
    $('molecule').dataset.frame = String(state.frame);
    $('molecule').dataset.path = state.source;
    $('selected-atoms').textContent = 'Highlighted: ' + (highlighted.map(atomLabel).join('–') || 'pick atoms in the molecule');
  }

  function svgElement(tag, attributes, text) {
    const node = document.createElementNS('http://www.w3.org/2000/svg',tag);
    Object.entries(attributes).forEach(([key,value]) => node.setAttribute(key,String(value)));
    if (text !== undefined) node.textContent = text;
    return node;
  }

  function renderChart(kind) {
    const entry = entryById(selected[kind]);
    const {svg, output} = chartParts[kind];
    svg.replaceChildren();
    svg.append(svgElement('title',{},entry.label + ': ' + kind + ' in ' + entry.unit + ', versus image index'));
    const width=520,height=148,left=53,right=13,top=13,bottom=31;
    const finite = [...entry.geodesic,...($('baseline').checked ? entry.cartesian : [])].filter(v => v !== null && Number.isFinite(v));
    let low=finite.length ? Math.min(...finite) : 0, high=finite.length ? Math.max(...finite) : 1;
    const pad = Math.max((high-low)*.13,kind==='distance' ? .015 : .8); low-=pad; high+=pad;
    const x = frame => left + frame / (count-1) * (width-left-right);
    const y = value => height-bottom-(value-low)/(high-low)*(height-top-bottom);
    for (let tick=0;tick<3;tick++) {
      const value=low+(high-low)*tick/2,py=y(value);
      svg.append(svgElement('line',{x1:left,y1:py,x2:width-right,y2:py,stroke:'#e4eaf1','stroke-width':1}));
      svg.append(svgElement('text',{x:left-9,y:py+4,'text-anchor':'end'},value.toFixed(kind==='distance'?2:0)));
    }
    [0,Math.round((count-1)/2),count-1].forEach(frame => {
      svg.append(svgElement('text',{x:x(frame),y:height-13,'text-anchor':'middle'},frame));
    });
    svg.append(svgElement('text',{x:width-right,y:height-1,'text-anchor':'end'},'Image index'));
    svg.append(svgElement('text',{x:4,y:11},entry.unit));
    const linePath = values => {
      let pen=false;
      return values.map((value,i) => {
        if(value===null || !Number.isFinite(value)){pen=false;return '';}
        const command=(pen?'L':'M')+x(i).toFixed(2)+','+y(value).toFixed(2);pen=true;return command;
      }).join(' ');
    };
    if ($('baseline').checked) svg.append(svgElement('path',{d:linePath(entry.cartesian),fill:'none',stroke:'#b56c25','stroke-width':1.7,'stroke-dasharray':'5 4'}));
    svg.append(svgElement('path',{d:linePath(entry.geodesic),fill:'none',stroke:'#087f8c','stroke-width':2.5}));
    svg.append(svgElement('line',{x1:x(state.frame),x2:x(state.frame),y1:top,y2:height-bottom,stroke:'#172944','stroke-width':1,'stroke-opacity':.5}));
    const current = entry[state.source][state.frame];
    if(current!==null && Number.isFinite(current)) svg.append(svgElement('circle',{
      cx:x(state.frame),cy:y(current),r:4,fill:state.source==='geodesic'?'#087f8c':'#b56c25',stroke:'white','stroke-width':1.5}));
    output.textContent = textValue(current,entry.unit);
    output.dataset.measurement = entry.id;
    svg.dataset.frame = state.frame;
    svg.onclick = event => {
      pause();
      const bounds=svg.getBoundingClientRect();
      const logical=(event.clientX-bounds.left)/bounds.width*width;
      state.frame=Math.max(0,Math.min(count-1,Math.round((logical-left)/(width-left-right)*(count-1))));
      state.active=entry.id; update();
    };
  }

  function refreshMenus() {
    Object.keys(names).forEach(kind => {
      const select=chartParts[kind].select;
      select.replaceChildren();
      data.measurements.filter(entry=>entry.kind===kind).forEach(entry=>{
        const option=document.createElement('option');option.value=entry.id;
        option.textContent=entry.label+(entry.id==='distance_1_4'?' · contact':'');
        select.append(option);
      });
      select.value=selected[kind];
    });
  }

  function update() {
    $('frame').value=state.frame;
    $('frame-value').textContent=state.frame;
    $('frame-caption').textContent='Image '+state.frame+' / '+(count-1);
    $('progress-value').textContent=Math.round(state.frame/(count-1)*100)+'% of image sequence';
    Object.keys(names).forEach(renderChart);
    renderMolecule();
  }

  function pause() {
    clearInterval(state.timer);state.timer=null;$('play').textContent='Play path';
  }

  function play() {
    if(state.timer){pause();return;}
    if(state.frame===count-1) state.frame=0;
    $('play').textContent='Pause';
    state.timer=setInterval(()=>{if(state.frame>=count-1){pause();return;}state.frame++;update();},1000/Number($('speed').value));
  }

  const sub=(a,b)=>a.map((v,i)=>v-b[i]);
  const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
  const norm=a=>Math.hypot(...a);
  const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
  function coordinate(xyz,atoms,kind) {
    if(kind==='distance')return norm(sub(xyz[atoms[0]],xyz[atoms[1]]));
    const a=sub(xyz[atoms[0]],xyz[atoms[1]]),b=sub(xyz[atoms[2]],xyz[atoms[1]]);
    if(kind==='angle') {
      const denominator=norm(a)*norm(b);
      return denominator>1e-14?Math.acos(Math.max(-1,Math.min(1,dot(a,b)/denominator)))*180/Math.PI:null;
    }
    const length=norm(b);if(length<=1e-14)return null;
    const axis=b.map(v=>v/length),c=sub(xyz[atoms[3]],xyz[atoms[2]]);
    const v=sub(a,axis.map(t=>t*dot(a,axis))),w=sub(c,axis.map(t=>t*dot(c,axis)));
    if(norm(v)<=1e-10*norm(a)||norm(w)<=1e-10*norm(c))return null;
    const value=Math.atan2(dot(cross(axis,v),w),dot(v,w))*180/Math.PI;
    return Math.abs(value+180)<1e-10?180:value;
  }
  function unwrap(values) {
    let previous=null,offset=0;
    return values.map(value=>{
      if(value===null){previous=null;offset=0;return null;}
      if(previous!==null){const delta=value-previous;if(delta>180)offset-=360;else if(delta< -180)offset+=360;}
      previous=value;return value+offset;
    });
  }

  $('custom-form').onsubmit=event=>{
    event.preventDefault();pause();
    const kind=$('custom-kind').value,required={distance:2,angle:3,dihedral:4}[kind];
    const tokens=$('custom-atoms').value.trim().split(/[\s,;]+/);
    const atoms=tokens.map(token=>Number(token)-1);
    if(tokens.some(token=>!/^\d+$/.test(token))||atoms.length!==required||new Set(atoms).size!==required||atoms.some(i=>i<0||i>=data.symbols.length)){
      $('custom-status').classList.add('error');$('custom-status').textContent='Choose '+required+' distinct atom numbers between 1 and '+data.symbols.length+'.';return;
    }
    const id=kind+'_'+atoms.map(i=>i+1).join('_');
    let entry=entryById(id);
    if(!entry){
      entry={id,kind,atoms,label:atoms.map(atomLabel).join('–'),unit:kind==='distance'?'Å':'°'};
      ['geodesic','cartesian'].forEach(source=>{
        const values=data.coordinates[source].map(xyz=>coordinate(xyz,atoms,kind));
        entry[source]=kind==='dihedral'?unwrap(values):values;
        if(kind==='dihedral')entry['wrapped_'+source]=values;
      });
      data.measurements.push(entry);
    }
    state.picking=false;state.picks=[];$('pick-atoms').textContent='Pick atoms';
    selected[kind]=id;state.active=id;refreshMenus();update();
    $('custom-status').classList.remove('error');$('custom-status').textContent='Tracking '+entry.label+' across '+count+' images. Included in the CSV download.';
  };

  function download(name,text,type='text/plain') {
    $('export-preview').hidden=false;
    $('export-name').textContent='Export text: '+name;
    $('export-text').value=text;
    const url=URL.createObjectURL(new Blob([text],{type}));
    const link=document.createElement('a');link.href=url;link.download=name;link.hidden=true;
    document.body.append(link);link.click();link.remove();
    setTimeout(()=>URL.revokeObjectURL(url),60000);
  }
  function xyz(frames) {
    return frames.map((frame,index)=>data.symbols.length+'\nimage '+index+'; Angstrom; '+state.source+'\n'+frame.map((point,i)=>data.symbols[i]+' '+point.map(v=>v.toFixed(10)).join(' ')).join('\n')).join('\n')+'\n';
  }
  $('download-path').onclick=()=>download('butane_'+state.source+'.xyz',xyz(data.coordinates[state.source]));
  $('download-frame').onclick=()=>download('butane_image_'+state.frame+'.xyz',xyz([data.coordinates[state.source][state.frame]]).replace('image 0;','image '+state.frame+';'));
  $('download-csv').onclick=()=>{
    const columns=data.measurements.flatMap(entry=>(entry.kind==='dihedral'?['geodesic','cartesian','wrapped_geodesic','wrapped_cartesian']:['geodesic','cartesian']).map(source=>({entry,source})));
    const header=['image',...columns.map(({entry,source})=>entry.id+'_'+source+(entry.unit==='Å'?'_A':'_deg'))];
    const rows=Array.from({length:count},(_,i)=>[i,...columns.map(({entry,source})=>entry[source][i]??'')].join(','));
    download('butane_measurements.csv',[header.join(','),...rows].join('\n')+'\n','text/csv');
  };

  Object.keys(names).forEach(kind=>{
    const block=document.createElement('section');block.className='chart-block';
    const heading=document.createElement('div');heading.className='chart-header';
    const label=document.createElement('label');label.textContent=names[kind];
    const select=document.createElement('select');select.id='select-'+kind;select.setAttribute('aria-label',names[kind]+' coordinate');label.append(select);
    const output=document.createElement('output');output.className='live-value';output.id='value-'+kind;
    const svg=svgElement('svg',{viewBox:'0 0 520 148',class:'geometry-chart',role:'img','aria-label':names[kind]+' along the interpolated path'});
    heading.append(label,output);block.append(heading,svg);$('charts').append(block);
    chartParts[kind]={select,svg,output};
    select.onchange=()=>{selected[kind]=select.value;state.active=select.value;update();};
  });
  $('frame').max=count-1;$('image-count').textContent=count;
  $('metadata').textContent=JSON.stringify(data.metadata,null,2);
  $('frame').oninput=()=>{pause();state.frame=Number($('frame').value);update();};
  $('first').onclick=()=>{pause();state.frame=0;update();};
  $('last').onclick=()=>{pause();state.frame=count-1;update();};
  $('play').onclick=play;
  $('speed').onchange=()=>{if(state.timer){pause();play();}};
  $('path-choice').onchange=()=>{state.source=$('path-choice').value;update();};
  $('baseline').onchange=()=>Object.keys(names).forEach(renderChart);
  $('show-h').onchange=()=>renderMolecule();$('show-labels').onchange=()=>renderMolecule();
  $('reset-view').onclick=()=>{if(viewer&&initialView){viewer.setView(initialView);viewer.render();}};
  $('pick-atoms').onclick=()=>{state.picking=!state.picking;state.picks=[];$('custom-atoms').value='';$('pick-atoms').textContent=state.picking?'Cancel picking':'Pick atoms';$('custom-status').textContent=state.picking?'Click atoms in the molecule in measurement order.':'Atom picking cancelled.';renderMolecule();};
  document.addEventListener('visibilitychange',()=>{if(document.hidden)pause();});
  try {
    viewer=$3Dmol.createViewer($('molecule'),{backgroundColor:'white',backgroundAlpha:0,antialias:true});
    new ResizeObserver(()=>{viewer.resize();viewer.render();}).observe($('molecule'));
  } catch(error) {
    $('viewer-error').hidden=false;$('viewer-error').textContent='The 3D viewer needs WebGL. Geometry graphs and downloads remain available.';
  }
  refreshMenus();update();
})();
