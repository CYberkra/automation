// Render the research Markdown with docx-js. No field data is read.
// node scripts/render_literature_review.cjs <absolute docx module directory>
const fs = require('fs');
const path = require('path');
const d = require(process.argv[2] || 'docx');
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  HeadingLevel, WidthType, AlignmentType, BorderStyle, ShadingType,
  Header, Footer, PageNumber, ExternalHyperlink, LevelFormat, TableLayoutType } = d;

const root = path.resolve(__dirname, '..');
const stem = '2026-09-23_landslide_literature_review';
const source = path.join(root, 'docs', 'research', stem + '.md');
const output = path.join(root, 'docs', 'research', stem + '.docx');
const font = { ascii: 'Calibri', hAnsi: 'Calibri', eastAsia: '宋体', cs: 'Calibri' };
const headingFont = { ascii: 'Calibri', hAnsi: 'Calibri', eastAsia: '微软雅黑', cs: 'Calibri' };
const ink = '17334A', accent = '176C84', muted = '607080', width = 9026;

function inline(s, props = {}) {
  const result = [];
  const pattern = /(\*\*([^*]+)\*\*|`([^`]+)`|\[([^\]]+)\]\(([^)]+)\)|\*([^*]+)\*)/g;
  let start = 0, m;
  while ((m = pattern.exec(s))) {
    if (m.index > start) result.push(new TextRun({text:s.slice(start,m.index), ...props}));
    if (m[2]) result.push(new TextRun({text:m[2], bold:true, ...props}));
    else if (m[3]) result.push(new TextRun({text:m[3], color:accent, ...props}));
    else if (m[4]) result.push(new ExternalHyperlink({link:m[5], children:[new TextRun({text:m[4], style:'Hyperlink', ...props})]}));
    else result.push(new TextRun({text:m[6], italics:true, ...props}));
    start = pattern.lastIndex;
  }
  if (start < s.length) result.push(new TextRun({text:s.slice(start), ...props}));
  return result;
}
function para(text, extra={}) { return new Paragraph({children:inline(text), keepLines:true, ...extra}); }
const border = {style:BorderStyle.SINGLE, size:4, color:'D8E2E8'};
function table(lines) {
  const values=lines.filter(x=>!/^\|\s*:?-+/.test(x)).map(x=>x.trim().slice(1,-1).split('|').map(t=>t.trim()));
  const n=values[0].length;
  const widths=n===4?(values[0][0]==='编号'?[760,2480,3360,2426]:[1510,2540,2730,2246]):n===3?[1700,3580,3746]:Array(n).fill(Math.floor(width/n));
  widths[n-1]=width-widths.slice(0,n-1).reduce((a,b)=>a+b,0);
  return new Table({width:{size:width,type:WidthType.DXA},columnWidths:widths,layout:TableLayoutType.FIXED,
    borders:{top:border,bottom:border,left:border,right:border,insideHorizontal:border,insideVertical:border},
    rows:values.map((row,i)=>new TableRow({tableHeader:i===0,cantSplit:true,children:row.map((cell,j)=>new TableCell({
      width:{size:widths[j],type:WidthType.DXA},
      margins:{top:95,bottom:95,left:105,right:105},
      shading:{type:ShadingType.CLEAR,fill:i===0?ink:(i%2?'F1F6F8':'FFFFFF')},
      children:[new Paragraph({children:inline(cell,{size:19,color:i===0?'FFFFFF':'243746',bold:i===0}),spacing:{after:0,line:285},widowControl:true})]
    }))}))});
}

const children = [
  new Paragraph({spacing:{before:900,after:300},children:[new TextRun({text:'UAV-GPR / SFCW',font:headingFont,size:24,color:accent,bold:true})]}),
  new Paragraph({spacing:{after:180},children:[new TextRun({text:'非显性滑坡场景下',font:headingFont,size:44,color:ink,bold:true})]}),
  new Paragraph({spacing:{after:360},children:[new TextRun({text:'自动处理文献综述',font:headingFont,size:48,color:ink,bold:true})]}),
  para('gprMax 仿真 · 背景抑制 · 增益补偿 · ML 配置选择',{spacing:{after:500},children:[new TextRun({text:'gprMax 仿真 · 背景抑制 · 增益补偿 · ML 配置选择',size:25,color:muted,font:headingFont})]}),
  para('研究综述 v1.0  |  2026 年 9 月 23 日',{spacing:{after:400}}),
  para('当前决定：场景特化、仿真优先、两类算子起步。',{spacing:{after:180}}),
  para('本轮阅读：8 篇论文正文；1 份会议摘要；2 份指南草案相关条款。',{spacing:{after:180}}),
  para('SGC 全文访问受限，证据缺口单列。所有方案尚待本项目验证。',{spacing:{after:500}}),
  para('阅读导航',{heading:HeadingLevel.HEADING_2}),
  para('第 1–4 节：结论、术语与逐篇证据'),
  para('第 5–6 节：仿真场景、参考定义与评价'),
  para('第 7–8 节：算子契约、ML 路线与待补证据'),
  para('第 9–10 节：原始来源与项目文档入口'),
  new Paragraph({pageBreakBefore:true,spacing:{after:0},children:[]}),
];

const lines=fs.readFileSync(source,'utf8').split(/\r?\n/);
let i=lines.findIndex(x=>x.startsWith('## 1.'));
let orderedInstance=0;
for(;i<lines.length;i++){
  const line=lines[i].trim(); if(!line)continue;
  if(line.startsWith('|')){const block=[];while(i<lines.length && lines[i].trim().startsWith('|')) block.push(lines[i++]);i--;if(block[0].includes('编号'))children.push(new Paragraph({pageBreakBefore:true,spacing:{after:0},children:[]}));children.push(table(block));children.push(new Paragraph({spacing:{after:40},children:[]}));continue;}
  const h=line.match(/^(#{2,3})\s+(.*)$/);
  if(h){const isMain=h[1].length===2;children.push(para(h[2],{heading:isMain?HeadingLevel.HEADING_1:HeadingLevel.HEADING_2}));continue;}
  const num=line.match(/^\d+\.\s+(.*)$/);
  if(num){if(!/^\d+\.\s/.test(lines[i-1]||''))orderedInstance++;children.push(para(num[1],{numbering:{reference:'ordered',level:0,instance:orderedInstance},spacing:{after:120}}));continue;}
  if(line.startsWith('- ')){children.push(para(line.slice(2),{numbering:{reference:'bullet',level:0},spacing:{after:120}}));continue;}
  const nextLine=lines.slice(i+1).find(x=>x.trim());
  children.push(para(line,{keepNext:!!nextLine && nextLine.trim().startsWith('|') && !nextLine.includes('编号')}));
}

const doc=new Document({
  creator:'UAV-GPR 自动处理研究',title:'非显性滑坡场景下 UAV-GPR 自动处理文献综述',description:'公开文献证据、仿真设计与两类处理评价；非性能报告',
  styles:{default:{document:{run:{font,size:21,color:'243746'},paragraph:{spacing:{after:140,line:330},widowControl:true}}},
    paragraphStyles:[
      {id:'Heading1',name:'heading 1',basedOn:'Normal',next:'Normal',quickFormat:true,run:{font:headingFont,size:32,bold:true,color:ink},paragraph:{spacing:{before:260,after:180},keepNext:true,keepLines:true,outlineLevel:0}},
      {id:'Heading2',name:'heading 2',basedOn:'Normal',next:'Normal',quickFormat:true,run:{font:headingFont,size:25,bold:true,color:accent},paragraph:{spacing:{before:220,after:140},keepNext:true,keepLines:true,outlineLevel:1}}
    ]},
  numbering:{config:[
    {reference:'bullet',levels:[{level:0,format:LevelFormat.BULLET,text:'\u2022',alignment:AlignmentType.LEFT,style:{paragraph:{indent:{left:280,hanging:220}}}}]},
    {reference:'ordered',levels:[{level:0,format:LevelFormat.DECIMAL,text:'%1.',alignment:AlignmentType.LEFT,style:{paragraph:{indent:{left:330,hanging:280}}}}]}
  ]},
  sections:[{properties:{page:{size:{width:11906,height:16838},margin:{top:1134,bottom:1134,left:1440,right:1440,header:550,footer:550}}},
    headers:{default:new Header({children:[new Paragraph({children:[new TextRun({text:'UAV-GPR 自动处理研究   /   文献综述',size:17,color:muted})],border:{bottom:{color:'C8D7DF',style:BorderStyle.SINGLE,size:4}},spacing:{after:100}})]})},
    footers:{default:new Footer({children:[new Paragraph({alignment:AlignmentType.RIGHT,children:[new TextRun({text:'研究方案 · 尚待验证    |    ',size:17,color:muted}),new TextRun({children:[PageNumber.CURRENT],size:17,color:muted}),new TextRun({text:' / ',size:17,color:muted}),new TextRun({children:[PageNumber.TOTAL_PAGES],size:17,color:muted})]})]})},children}]
});
Packer.toBuffer(doc).then(buffer=>{fs.writeFileSync(output,buffer);process.stdout.write(output+'\n');});
