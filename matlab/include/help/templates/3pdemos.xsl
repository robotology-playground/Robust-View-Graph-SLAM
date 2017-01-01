<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet
  version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:file="java.io.File">

<!--
This is an XSL stylesheet which converts DEMOS.XML files into HTML.

Copyright 1984-2015 The MathWorks, Inc.
-->

  <xsl:output method="html" indent="yes"/>

  <xsl:param name="demosroot"/>
  <xsl:param name="languageDir"/>
  <xsl:param name="mfile"/>
  <xsl:param name="mfiledesc"/>
  <xsl:param name="mgui"/>
  <xsl:param name="mguidesc"/>
  <xsl:param name="model"/>
  <xsl:param name="modeldesc"/>
  <xsl:param name="productlink"/>
  <xsl:param name="section"/>
  <xsl:param name="subsection"/>
  <xsl:param name="uses"/>
  <xsl:param name="video"/>
  <xsl:param name="videodesc"/>
  <xsl:param name="matlabroot"/>
  <xsl:param name="docroot"/>

<!--xsl:variable name="filesep" select="system-property('file.separator')"/-->
<xsl:variable name="filesep" select="'/'"/><!-- TODO: Pass this system property from Java? -->
<xsl:variable name="language" select="$languageDir"/>

<!-- Windows requires three slashes and Unix doesn't mind them. -->
<xsl:variable name="urlbase">file:///</xsl:variable>

<xsl:variable name="demosrooturl">
  <xsl:value-of select="$urlbase"/>
  <xsl:call-template name="globalReplace">
    <xsl:with-param name="outputString" select="$demosroot"/>
    <xsl:with-param name="target" select="'\'"/>
    <xsl:with-param name="replacement" select="'/'"/>
  </xsl:call-template>
</xsl:variable>

<!-- Create the full product name. -->
<xsl:variable name="productName">
  <xsl:value-of select="/demos/name"/>
  <xsl:choose>
    <xsl:when test="/demos/type='toolbox'"> Toolbox</xsl:when>
    <xsl:when test="/demos/type='blockset'"> Blockset</xsl:when>
  </xsl:choose>
</xsl:variable>

<xsl:template match="demos">

<!--html xmlns="http://www.w3.org/1999/xhtml">
<head>
  <title><xsl:value-of select="$productName"/> Examples</title>
  <base href="{$demosrooturl}/"/>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
  <link rel="stylesheet" type="text/css"
    href="{$urlbase}{$matlabroot}/toolbox/matlab/helptools/private/demoindex.css"/>
</head>

<body-->

<div id="background-container">
<!-- BODY BEGIN -->


    <!-- banner -->
  <div id="header-container">
    <div id="header">
      <div>
        <span class="header-demo-name"><xsl:value-of select="$productName"/></span>&#160;
        <span class="header-demos">EXAMPLES</span>
      </div>

    <!-- description -->
    <div class="header-content">
      <xsl:choose>
        <xsl:when test="description[@isCdata='no']">
          <xsl:copy-of select="description/*"/>
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="description" disable-output-escaping="yes"/>
        </xsl:otherwise>
      </xsl:choose>
    </div>

      <!-- link to website -->
      <div class="header-link">
        <xsl:choose>

          <!-- Allow non-MathWorks products to override with their own link specified in demos.xml. -->
          <xsl:when test="website">
            <xsl:copy-of select="website/*"/>
          </xsl:when>
          <xsl:otherwise>

            <!-- Determine the product's unique directory to act as an identifier. -->

            <xsl:variable name="pathAfterToolbox">
              <xsl:value-of select="substring-after($demosroot,concat('toolbox',$filesep))"/>
            </xsl:variable>

            <xsl:variable name="startsWithUniqueDir">
               <xsl:choose>
                 <xsl:when test="$pathAfterToolbox=concat('simulink',$filesep,'simdemos')">
                   <!-- "simulink/simdemos" to "simulink" (special case of next rule) -->
                   <xsl:value-of select="'simulink'"/>
                 </xsl:when>
                 <xsl:when test="starts-with($pathAfterToolbox,'simulink') or starts-with($pathAfterToolbox,'physmod')">
                   <!-- "simulink/fixedandfloat" to "fixedandfloat" (down one) -->
                   <xsl:value-of select="substring-after($pathAfterToolbox,$filesep)"/>
                 </xsl:when>
                 <xsl:when test="$pathAfterToolbox=concat('rtw',$filesep,'rtwdemos')">
                   <!-- "rtw/rtwdemos" to "rtw" (special case of next rule) -->
                   <xsl:value-of select="'rtw'"/>
                 </xsl:when>
                 <xsl:when test="starts-with($pathAfterToolbox,'rtw')">
                   <!-- "rtw/targets/hc12/hc12" to "hc12" (down two) -->
                   <xsl:value-of select="substring-after(substring-after($pathAfterToolbox,$filesep),$filesep)"/>
                 </xsl:when>
                 <xsl:otherwise>
                   <xsl:value-of select="$pathAfterToolbox"/>
                 </xsl:otherwise>
               </xsl:choose>
            </xsl:variable>

            <xsl:variable name="uniqueProductDir">
              <xsl:choose>
                <xsl:when test="contains($startsWithUniqueDir,$filesep)">
                  <!-- "bioinfo/biodemos" to "bioinfo" -->
                  <xsl:value-of select="substring-before($startsWithUniqueDir,$filesep)"/>
                </xsl:when>
                <xsl:otherwise>
                  <!-- "stats" to "stats" -->
                  <xsl:value-of select="$startsWithUniqueDir"/>
                </xsl:otherwise>
              </xsl:choose>
            </xsl:variable>

            <!-- Add the link to the page. -->
            <a href="http://www.mathworks.com/pl_{$uniqueProductDir}"><xsl:value-of select="$productlink"/></a>&#160;
            <a href="http://www.mathworks.com/pl_{$uniqueProductDir}">
                <img src="{$urlbase}{$matlabroot}/toolbox/matlab/helptools/private/demo_arrow.png" border="none" alt="arrow" align="middle" />
            </a>
          </xsl:otherwise>
        </xsl:choose>
      </div>

    </div>
  </div>

<!-- content -->
<div class="demopanel-container">
<div class="demopanel">
<xsl:choose>
  <xsl:when test="not($section)">
    <!-- Show all sections (if any). -->
    <xsl:apply-templates select="demoitem | demosection"/>
  </xsl:when>
  <xsl:otherwise>
    <!-- Show only one category. -->
    <xsl:choose>
      <xsl:when test="$subsection">
        <xsl:apply-templates select="//demosection[concat(translate(label,' ','-'),count(preceding::demosection))=$section]/demosection[label=$subsection]"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:apply-templates select="//demosection[concat(translate(label,' ','-'),count(preceding::demosection))=$section]"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:otherwise>
</xsl:choose>
</div>
</div>
</div>
<!--/body>
</html-->

</xsl:template>


<!-- Display a demosection. -->
<xsl:template match="demosection">

  <!-- section heading -->
    <xsl:choose>
      <xsl:when test="count(ancestor::*)=1">
        <!-- Section -->
        <div class="demopanel">
          <div class="demopanel-title">
            <xsl:value-of select="label"/>
          </div>
          <xsl:apply-templates select="demoitem | demosection"/>
        </div>
      </xsl:when>
      <xsl:when test="count(ancestor::*)=2">
        <!-- Sub-section -->
        <div>
          <xsl:attribute name="class">
            <xsl:choose>
              <xsl:when test="../demoitem">
                demopanel-sub-container
              </xsl:when>
              <xsl:otherwise>
                demopanel-container
              </xsl:otherwise>
            </xsl:choose>
          </xsl:attribute>
          <div class="demopanel-sub-title">
            <xsl:value-of select="label"/>
          </div>
          <xsl:apply-templates select="demoitem | demosection"/>
        </div>
      </xsl:when>
    </xsl:choose>

</xsl:template>


<!-- Display a demo. -->
<xsl:template match="demoitem">

  <xsl:variable name="link">
    <!-- Construct title with link. -->
    <xsl:choose>
      <xsl:when test="type and file">
		<xsl:variable name="localFile">
			<xsl:choose>
				<xsl:when test="$language != ''">
					<xsl:call-template name="buildLocalizedPath"><xsl:with-param name="filepath" select="file" /></xsl:call-template>
				</xsl:when>
				<xsl:otherwise></xsl:otherwise>
			</xsl:choose>
		</xsl:variable>						
		<xsl:call-template name="checkLocalizedFile">
			<xsl:with-param name="localizedFile" select="$localFile" />
			<xsl:with-param name="defaultFile" select="file" />
		</xsl:call-template>
      </xsl:when>
      <xsl:when test="'M-file'=type">
        <xsl:variable name="source"> <xsl:value-of select="source"/> </xsl:variable>
		<xsl:variable name="localFile">
			<xsl:choose>
				<xsl:when test="$language != ''"><xsl:value-of select="concat('html/',$language,'/',$source,'.html')" /></xsl:when>
				<xsl:otherwise></xsl:otherwise>
			</xsl:choose>
		</xsl:variable>						
        <xsl:variable name="defaultFile" select="concat('html/',$source,'.html')" />
        <xsl:call-template name="checkLocalizedFile">
          <xsl:with-param name="localizedFile" select="$localFile" />
          <xsl:with-param name="defaultFile" select="$defaultFile" />
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:text>matlab:indexhelper('</xsl:text>
        <xsl:call-template name="quotify">
          <!-- Strip off the leading file:/// for use with indexhelper. -->
<xsl:with-param name="string" select="substring-after($demosroot,'file:///')"/>
        </xsl:call-template>
        <xsl:text>','</xsl:text>
        <xsl:call-template name="quotify">
          <xsl:with-param name="string" select="source"/>
        </xsl:call-template>
        <xsl:text>','</xsl:text>
        <xsl:call-template name="quotify">
          <xsl:with-param name="string" select="callback"/>
        </xsl:call-template>
        <xsl:text>','</xsl:text>
        <xsl:call-template name="quotify">
          <xsl:with-param name="string" select="$productName"/>
        </xsl:call-template>
        <xsl:text>','</xsl:text>
        <xsl:call-template name="quotify">
          <xsl:with-param name="string" select="label"/>
        </xsl:call-template>
        <xsl:text>','</xsl:text>
        <xsl:call-template name="quotify">
          <xsl:with-param name="string" select="file"/>
        </xsl:call-template>
        <xsl:text>')</xsl:text>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:variable>

  <!-- Create the HTML for one demoitem -->
  <table width="100%" border="0" cellspacing="0" cellpadding="0">
    <tr>
      <!-- thumbnail -->
      <!-- TODO: Remove for now since we cannot check for file existence -->
      <!--td class="demopanel-thumbnail">
        <a href="{$link}">
          <img>
            <xsl:attribute name="src">
              <xsl:call-template name="globalReplace">
                <xsl:with-param name="outputString">
                  <xsl:call-template name="chooseThumbnail"/>
                </xsl:with-param>
                <xsl:with-param name="target" select="'\'"/>
                <xsl:with-param name="replacement" select="'/'"/>
              </xsl:call-template>
            </xsl:attribute>
          </img>
        </a>
      </td-->

      <!-- text -->
      <td class="demopanel-description">
        <xsl:attribute name="colspan">
          <xsl:choose>
            <xsl:when test="count(ancestor::*)=1">2</xsl:when>
            <xsl:when test="count(ancestor::*)=2">2</xsl:when>
            <xsl:when test="count(ancestor::*)=3">1</xsl:when>
          </xsl:choose>
        </xsl:attribute>

        <xsl:variable name="linkPath">
          <xsl:call-template name="checkLinkPath">
            <xsl:with-param name="linkPath" select="$link"/>
            <xsl:with-param name="rootPath" select="$demosroot"/>
          </xsl:call-template>
        </xsl:variable> 
        <a href="{$linkPath}"><xsl:value-of select="label"/></a>
        <br/>

        <xsl:variable
          name="productDependencies"
          select="dependency[(.!='student') and (.!='nonstudent')]"/>
        <xsl:if test="$productDependencies">
          <span class="dependency"><xsl:value-of select="$uses"/>&#160;<xsl:apply-templates select="$productDependencies"/></span>
        </xsl:if>

      </td>

      <!-- type -->
      <td width="80" class="demopanel-type" align="left">
        <xsl:choose>
          <xsl:when test="type='M-file'">
            <xsl:attribute name="title"><xsl:value-of select="$mfiledesc"/></xsl:attribute>
            <img src="{$urlbase}{$matlabroot}/toolbox/matlab/helptools/private/type_m-file.png" align="middle"/><xsl:value-of select="$mfile"/>
          </xsl:when>
          <xsl:when test="type='M-GUI'">
            <xsl:attribute name="title"><xsl:value-of select="$mguidesc"/></xsl:attribute>
            <img src="{$urlbase}{$matlabroot}/toolbox/matlab/helptools/private/type_m-gui.png" align="middle"/><xsl:value-of select="$mgui"/>
          </xsl:when>
          <xsl:when test="type='model'">
            <xsl:attribute name="title"><xsl:value-of select="$modeldesc"/></xsl:attribute>
            <img src="{$urlbase}{$matlabroot}/toolbox/matlab/helptools/private/type_model.png" align="middle"/><xsl:value-of select="$model"/>
          </xsl:when>
          <xsl:when test="type='video'">
            <xsl:attribute name="title"><xsl:value-of select="$videodesc"/></xsl:attribute>
            <img src="{$urlbase}{$matlabroot}/toolbox/matlab/helptools/private/type_video.png" align="middle"/><xsl:value-of select="$video"/>
          </xsl:when>
        </xsl:choose>
      </td>
    </tr>
  </table>

</xsl:template>


<xsl:template match="dependency">
  <xsl:if test="position()!=1">, </xsl:if>
  <xsl:value-of select="."/>
</xsl:template>

<!-- Figure out which thumbnail to use for this demo. -->
<xsl:template name="chooseThumbnail">
  <xsl:variable name="thumbnail" select="concat($demosroot,'/',substring-before(file,'.html'),'_img_thumbnail.png')"/>
  <xsl:variable name="thumbnail2" select="concat($demosroot,'/html/',callback,'.png')"/>
  <xsl:variable name="thumbnail3" select="concat($demosroot,'/html/',m-file,'.png')"/>
  <xsl:variable name="thumbnail4" select="concat($demosroot,'/',substring-before(file,'.html'),'.png')"/>
  <xsl:variable name="thumbnail5" select="concat($demosroot,'/html/',source,'.png')"/>
  <xsl:value-of select="$urlbase"/>
  <xsl:choose>
    <xsl:when test="file:exists(file:new($thumbnail5))"><xsl:value-of select="$thumbnail5"/></xsl:when>
    <xsl:when test="file:exists(file:new($thumbnail))"><xsl:value-of select="$thumbnail"/></xsl:when>
    <xsl:when test="file:exists(file:new($thumbnail2))"><xsl:value-of select="$thumbnail2"/></xsl:when>
    <xsl:when test="file:exists(file:new($thumbnail3))"><xsl:value-of select="$thumbnail3"/></xsl:when>
    <xsl:when test="file:exists(file:new($thumbnail4))"><xsl:value-of select="$thumbnail4"/></xsl:when>
    <xsl:when test="type='M-file'"><xsl:value-of select="$matlabroot"/>/toolbox/matlab/helptools/private/type_m-file_lg.png</xsl:when>
    <xsl:when test="type='video'"><xsl:value-of select="$matlabroot"/>/toolbox/matlab/helptools/private/type_video_lg.png</xsl:when>
    <xsl:when test="/demos/icon">
      <xsl:choose>
        <xsl:when test="starts-with(/demos/icon,'$matlabroot')"><xsl:value-of select="concat($matlabroot,substring-after(/demos/icon,'$matlabroot'))"/></xsl:when>
        <xsl:when test="starts-with(/demos/icon,'$toolbox')"><xsl:value-of select="concat($matlabroot,'/toolbox',substring-after(/demos/icon,'$toolbox'))"/></xsl:when>
        <xsl:otherwise><xsl:value-of select="/demos/icon"/></xsl:otherwise>
      </xsl:choose>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="'{$matlabroot}/toolbox/matlab/icons/matlabicon.gif'"/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>


<xsl:template name="quotify">
  <xsl:param name="string"/>
  <xsl:variable name="percentEscaped">
    <xsl:call-template name="globalReplace">
      <xsl:with-param name="outputString" select="$string"/>
      <xsl:with-param name="target" select="'%'"/>
      <xsl:with-param name="replacement" select="'%25'"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="doubleUnQuoted">
    <xsl:call-template name="globalReplace">
      <xsl:with-param name="outputString" select="$percentEscaped"/>
      <xsl:with-param name="target" select="'&quot;'"/>
      <xsl:with-param name="replacement" select="'&amp;quot;'"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:call-template name="globalReplace">
    <xsl:with-param name="outputString" select="$doubleUnQuoted"/>
    <xsl:with-param name="target" select='"&apos;"'/>
    <xsl:with-param name="replacement" select='"&apos;&apos;"'/>
  </xsl:call-template>
</xsl:template>


<!-- Search and replace  -->
<!-- From http://www.xml.com/lpt/a/2002/06/05/transforming.html -->
<xsl:template name="globalReplace">
  <xsl:param name="outputString"/>
  <xsl:param name="target"/>
  <xsl:param name="replacement"/>
  <xsl:choose>
    <xsl:when test="contains($outputString,$target)">
      <xsl:value-of select=
      "concat(substring-before($outputString,$target),$replacement)"/>
      <xsl:call-template name="globalReplace">
        <xsl:with-param name="outputString"
          select="substring-after($outputString,$target)"/>
        <xsl:with-param name="target" select="$target"/>
        <xsl:with-param name="replacement"
          select="$replacement"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="$outputString"/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<!-- A helper method used to build a localized path from a non-localized path. 
	 This method is called only when $language is not '' and we know that a 
	 localized directory should exist. -->
<xsl:template name="buildLocalizedPath">
   <xsl:param name="filepath" />
   <xsl:choose>
      <xsl:when test="contains($filepath, '/')"><xsl:value-of
        select="concat(substring-before($filepath,'/'), '/')" /><xsl:call-template name="buildLocalizedPath">
        <xsl:with-param name="filepath" select="substring-after($filepath,'/')" />
        </xsl:call-template></xsl:when>
      <xsl:otherwise>
		<xsl:value-of select="$language" />/<xsl:value-of select="$filepath" />
	  </xsl:otherwise>
   </xsl:choose>
</xsl:template>

<!-- A helper function that checks if a localized file path exists.  If it
     does, that file path will be returned, otherwise a specified default
     path will be. -->
<xsl:template name="checkLocalizedFile">
    <xsl:param name="localizedFile" />
    <xsl:param name="defaultFile" />
	<xsl:choose>
		<xsl:when test="$localizedFile != ''">
			<xsl:variable name="fullpath" select="concat($demosroot,'/',$localizedFile)" />
			<xsl:choose>
				<xsl:when test="file:exists(file:new($fullpath))"><xsl:value-of select="$localizedFile" /></xsl:when>
				<xsl:otherwise><xsl:value-of select="$defaultFile" /></xsl:otherwise>
			</xsl:choose>		
		</xsl:when>
		<xsl:otherwise><xsl:value-of select="$defaultFile" /></xsl:otherwise>
	</xsl:choose>	
</xsl:template>

<!-- A helper function to check the link path. If the link path is not a 
     matlab: link, pre-prend the root path. -->
<xsl:template name="checkLinkPath">
	<xsl:param name="linkPath"/>    
  <xsl:param name="rootPath"/>
	<xsl:choose>
		<xsl:when test="contains($linkPath,'matlab:')"><xsl:value-of select="$linkPath" /></xsl:when>      
		<xsl:otherwise><xsl:value-of select="concat($rootPath,'/',$linkPath)"/></xsl:otherwise>      
	</xsl:choose>    
</xsl:template>

</xsl:stylesheet>
